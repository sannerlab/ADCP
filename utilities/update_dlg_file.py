#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 14:28:57 2023

@author: sshanker
To make old dlg file format homegeneous with latest format. 
"""

def get_docking_data_from_dlg(dlgfilein):
    fid = open(dlgfilein,'r')
    dlgdata = fid.readlines()    
    line_read = False
    ignore_lines = 3
    header_lines = []
    score_lines =[]
    for l in dlgdata:
        if 'mode |  affinity' in l:
            line_read = True
        if line_read:
            if ignore_lines > 0:
                header_lines.append(l.rstrip())
                ignore_lines -= 1
            else:
                if not l.split()[0].isnumeric():
                    line_read = False
                    continue
                score_lines.append(l.rstrip())
   
    if len(score_lines) < 1: # in case no scores are there 
        return ['',]*3, None   
    return header_lines,score_lines, dlgdata

import numpy as np
def calculate_reranking_index(enzs_array, rerank_by_Eint=False, rearrangeposes = False):  
    #for rearraging poses
    rerank_metric = []  
    ommrank_adcprank_data = []
    rerank_method =""
    # [0,flnm[:-4]+"_fixed_min.pdb", enzs, not_restrained_res]
    for i in enzs_array:
        if rerank_by_Eint: # Ecomp-Erec-Epep
            rerank_metric.append(i[1][-2])    
                  
        else:
            rerank_metric.append(i[1][-1])     
    reaarange_for_omm_ranking = np.argsort(rerank_metric)  
    if rearrangeposes == False:        
              
        for current_v, rearrage_idx in enumerate(reaarange_for_omm_ranking):
            #modelnumber, newranking, adcp ranking, scores
            ommrank_adcprank_data.append([ current_v+1, current_v+1, 
                                          enzs_array[rearrage_idx]])
        #ommrank_adcprank_data = np.array(ommrank_adcprank_data)
        adcp_ranks = np.array([ i[1] for i in ommrank_adcprank_data])   

    else:
        reaarange_for_omm_ranking= list(reaarange_for_omm_ranking)
        for i, enz in enumerate(enzs_array):
            ommrank_adcprank_data.append([i+1,reaarange_for_omm_ranking.index(i)+1, enz])
        # self.myprint("\nOMM Ranking:NOT REARRANGING output poses using OpenMM (%s) energy" % rerank_method)
        # self.rearranged_data_as_per_asked = [ ommrank_adcprank_data[i] for i in adcp_ranks.argsort()]
    return ommrank_adcprank_data


def get_omm_data_from_dlgdata(summary_data, rerank_by_Eint=False, rearrangeposes = False):       
    omm_data = []
    enzs_array = []
    for i, line in enumerate(summary_data):
        if line.startswith('OMM Energy: Working on'):
            model_num = int(line.split("#")[1].split()[0])
            enzs=[ float(v) for v in summary_data[i+1].replace("="," ").replace(";","").split()[1:] if not v.startswith('E')]
            enzs.append(enzs[0] - enzs[1] -enzs[2])   # interaction energy
            enzs.append(enzs[0] - enzs[1])            # complex - protein energy
            enzs = [model_num, enzs]                  
            enzs_array.append(enzs)               
    
    omm_data = calculate_reranking_index(enzs_array, rerank_by_Eint=rerank_by_Eint, rearrangeposes = rearrangeposes) # [modelnumer, rankomm, [rankadcp, enzs]]        
    return omm_data
                
def update_dlg(**kw):  
    dlgfilein = kw['dlgfile']
    rerank_by_Eint = kw['reint']
    rearrangeposes = kw['dockingRanking']
    
    docking_header, docking_score, all_dlg_data = get_docking_data_from_dlg(dlgfilein)    
    omm_data = get_omm_data_from_dlgdata(all_dlg_data,rerank_by_Eint,rearrangeposes)
    if omm_data == None:        
        print("Already the lastest format")
        
        return
    if not len(omm_data):
        print("No previous minimization data is found. Docking data format is up-to-date.")
        # return
        
    new_data = []
    if rerank_by_Eint:
        rerank_method = "Ecomplex -Ereceptor -Epeptide"  
    else:
        rerank_method = "Ecomplex -Ereceptor"
    for line in all_dlg_data:
        if line.startswith('OMM Ranking:'):
            if "Ecomplex" in line:
                rerank_method = line.split("(")[1].split(")")[0]
            break                     
        new_data.append(line)
    pdmc=""
    
    for line in all_dlg_data[-3:]:
        if line.startswith("Post Docking Minimization Command"):
            pdmc = line
            break           
    if rearrangeposes:
        new_data.append('OMM Ranking:NOT REARRANGING output poses using OpenMM (%s) energy\n' % rerank_method)
    else:
        new_data.append("OMM Ranking:REARRANGING output poses using OpenMM (%s) energy\n" % rerank_method)   
    docking_data_index_after_mode = docking_header[2].find("+")+1
    new_data.append("OMM Ranking:                     +<-------OMMscore-------->+<-----------AutoDock CrankPep Scores"+
                 "-"*(len(docking_header[2])-docking_data_index_after_mode-38)+">+\n")
    new_data.append("OMM Ranking:-------+------+------+------------+------------+" + docking_header[2][docking_data_index_after_mode:]+"\n")        
    new_data.append("OMM Ranking: Model | Rank | Rank | E_Complex  |  E_Complex |"  + docking_header[0][docking_data_index_after_mode:]+"\n")
    new_data.append("OMM Ranking: #     |OpenMM| ADCP |-E_Receptor |-E_rec-E_pep|" + docking_header[1][docking_data_index_after_mode:]+"\n")
    new_data.append("OMM Ranking:-------+------+------+------------+------------+" + docking_header[2][docking_data_index_after_mode:]+"\n") 
    
    for model_num, omm_rank, adcpNenz in omm_data:
        adcp_rank, enz = adcpNenz
        ecomp,erec,epep,eint,ecomres = enz
        docking_score_line = docking_score[adcp_rank-1][docking_data_index_after_mode:]        
        new_data.append("OMM Ranking: %6d %6d %6d %10.1f   %10.1f   " %(model_num , omm_rank, adcp_rank, ecomres, eint ) + docking_score_line+"\n")                    
    new_data.append("OMM Ranking:-------+------+------+------------+------------+"+ docking_header[2][docking_data_index_after_mode:]+"\n")

    fid = open(dlgfilein+".dlgu",'w+')
    for nline in new_data:
        fid.write(nline)
    if pdmc.__len__():
        fid.write(pdmc.strip())
    fid.close()
    print("Created new file %s with updated format!" % (dlgfilein+".dlgu"))
    print("Check it and rename it if required!")

    
if __name__=='__main__':
    import sys
    import argparse
    parser = argparse.ArgumentParser(description='ADCP dlg file manipulator', 
                                  usage="usage: python %(prog)s -i dlgfile")
                                  # version="%prog 0.1")
    parser.add_argument('--version', action='version', version="0.0.1" )
    parser.add_argument("-i", "--dlgfile",dest="dlgfile",
                        help="input dlg file for manipulation")    

    parser.add_argument("-dr", "--dockingRanking", action="store_true",default=False,
                       dest="dockingRanking", help=( 'When docking solutions are\
                       minimized, the -dr flag is used to prevent the re-ordering\
                       and report the minimized solutions based on docking score\
                       ordered from best to worse. Default is False.'))                            
                             
    parser.add_argument("-reint", "--rerankByEinteraction",dest="reint",                         
                        action="store_true", default=False,
                        help=("To rank poses by openMM interaction energy (Ecomplex -Ereceptor -Epeptide). \
                              By default pose ranking will be performed by Ecomplex -Ereceptor."))  
    
                                    #OMM new line
    if len(sys.argv) == 1:
        parser.print_help()
        
    else:   
        kwn, unknw = parser.parse_known_args()
        if len(unknw)> 0:
            print('Unknown options: ',end="")
            for u in unknw:
                if u.startswith("-"):
                    print(u, end=" ")
        else:
            kw = vars(kwn)
            update_dlg(**kw)

    
    
    
