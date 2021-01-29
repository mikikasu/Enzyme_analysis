# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 23:08:47 2020

@author: kasug
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from natsort import natsorted
import collections
import os

def makeSeqRecord(seq_obj, id_for_fasta):
    record = SeqRecord(seq_obj, id=id_for_fasta, description="")
    return record

def make_fasta(gb_file,out_fasta,count):
    FASTA_id_withEC = []
    CDS_withEC_count = 0
    CDS_withaEC_count = 0
    CDS_withmultipleEC_count = 0
    
    
    for record in SeqIO.parse(gb_file, 'genbank'):
        # basket to collect CDS slices
        name = record.annotations["organism"]
        pos = name.find(' ')
        species = name[pos+1:]
        species_name = name[0] + '_' + species
        species_name = species_name[:6]
        
        CDS_number = 0
        CDS_withoutproduct_number = 0
    
    
        for feature in record.features:
            # pick only CDS
            if feature.type == 'CDS':
                # get EC number annotation
                CDS_number += 1
                
                if 'EC_number' in feature.qualifiers:
                    CDS_withEC_count += 1
                    #get EC_number's list
                    EC_list = feature.qualifiers['EC_number']   
                    
                    # To Extruct features whose at least one EC number includes hyphen
                    EC_list_in =[s for s in EC_list if '-' in s]
                    if len(EC_list_in) == 0 :
                        
                        if 'protein_id' in feature.qualifiers:
                            if 'gene' in feature.qualifiers:
                               gene = feature.qualifiers['gene'] 
                               gene_str = ''.join(gene)
                            else:
                                gene_str = 'No_gene'
                            
                            protein_id = feature.qualifiers['protein_id'] 
                            product = feature.qualifiers['product'] 
                           
                            EC_list = natsorted(EC_list)

                            if len(EC_list) >= 2 :
                                EC_list_3 = [ i[:5] for i in EC_list]
                                EC_list_3_unique = list(dict.fromkeys(EC_list_3))
                                c = collections.Counter(EC_list_3)
                                EC_list_3_multiple = [i[0] for i in c.items() if i[1] >= 2]
                            
                                if len(EC_list) != len(EC_list_3_unique):
                                    new_EC_list = EC_list
                                    for ECnumber_3 in EC_list_3_multiple:
                                        new_EC_list = [s for s in new_EC_list if (ECnumber_3 in s) == False]
                                        ECnumber_3_multiple = ECnumber_3 + '.X'
                                        new_EC_list.append(ECnumber_3_multiple)
                                        new_EC_list = natsorted(new_EC_list)
                                
                                    if len(new_EC_list) > 3:
                                        new_EC_list = new_EC_list[:3]
                                        EC_list_str = '+' + '_'.join(new_EC_list) 
                                    else:
                                        EC_list_str = '_'.join(new_EC_list)

                                else:
                                    if len(EC_list) < 3:
                                        EC_list_str = '_'.join(EC_list)
                                    else:
                                        new_EC_list = EC_list[:3]
                                        EC_list_str = '+' + '_'.join(new_EC_list)
                            else:
                                EC_list_str = '_'.join(EC_list)
                            
                                                            
                            EC_list_discription = ' '.join(EC_list)
                            protein_id_str = '_'.join(protein_id)
                            product_str = ''.join(product)  
 

                            if  'translation' in feature.qualifiers:
                                translation =  feature.qualifiers['translation'] 
                            else:
                                continue
                            
                            if len(EC_list_str) < 10:
                                #FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name +'_'+ protein_id_str +'_' + 'CDS' + str(CDS_number)+ ' ' + EC_list_discription +' '+ gene_str  + ' ' + 'CDS' + str(CDS_number) + '\n'+translation[0] + '\n'        
                                FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name +'_'+ protein_id_str +'_' + 'CDS' + str(CDS_number)+ '\n'+translation[0] + '\n'        
                                
                                #FASTA_id_withEC_str ='>'+ EC_list_str+'_'+ gene_str +'_'+ product_str + ' ' + protein_id_str + ' ' + 'CDS' + str(CDS_number)+'\n'+translation[0]          
                            else:
                                #FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name +'_'+ protein_id_str + ' ' + EC_list_discription +' '+ gene_str  + ' ' + 'CDS' + str(CDS_number) + '\n'+translation[0] + '\n'        
                                FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name +'_'+ protein_id_str + '\n'+translation[0] + '\n'        
                                
                                #FASTA_id_withEC_str ='>'+ EC_list_str+'_'+ gene_str +'_'+ product_str + ' ' + protein_id_str + ' ' + 'CDS' + str(CDS_number)+'\n'+translation[0]          
                            
                            
                            with open(out_fasta, 'a') as f:
                                print(FASTA_id_withEC_str, file = f)
                            FASTA_id_withEC.append(FASTA_id_withEC_str)
    
                        else:
                            if 'gene' in feature.qualifiers:
                               gene = feature.qualifiers['gene'] 
                               gene_str = ''.join(gene)
                            else:
                                gene_str = 'No_gene'
                                
                            protein_id = [None]
                            product = feature.qualifiers['product'] 
                            
                            EC_list = natsorted(EC_list)

                            if len(EC_list) >= 2 :
                                EC_list_3 = [ i[:5] for i in EC_list]
                                EC_list_3_unique = list(dict.fromkeys(EC_list_3))
                                c = collections.Counter(EC_list_3)
                                EC_list_3_multiple = [i[0] for i in c.items() if i[1] >= 2]
                            
                                if len(EC_list) != len(EC_list_3_unique):
                                    new_EC_list = EC_list
                                    for ECnumber_3 in EC_list_3_multiple:
                                        new_EC_list = [s for s in new_EC_list if (ECnumber_3 in s) == False]
                                        ECnumber_3_multiple = ECnumber_3 + '.X'
                                        new_EC_list.append(ECnumber_3_multiple)
                                        new_EC_list = natsorted(new_EC_list)
                                
                                    if len(new_EC_list) > 3:
                                        new_EC_list = new_EC_list[:3]
                                        EC_list_str = '+' + '_'.join(new_EC_list) 
                                    else:
                                        EC_list_str = '_'.join(new_EC_list)

                                else:
                                    if len(EC_list) < 3:
                                        EC_list_str = '_'.join(EC_list)
                                    else:
                                        new_EC_list = EC_list[:3]
                                        EC_list_str = '+' + '_'.join(new_EC_list)
                            else:
                                EC_list_str = '_'.join(EC_list)

                                                            
                            EC_list_discription = ' '.join(EC_list)                                
                            product_str = ''.join(product) 
                            
                            if  'translation' in feature.qualifiers:
                                translation =  feature.qualifiers['translation'] 
                            else:
                                continue
                            if len(EC_list_str) < 10:                            
                                #FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name +'_'+'NO protein_id'+'_' + 'CDS' + str(CDS_number) +' '+ gene_str  + ' ' + EC_list_discription + ' ' + 'CDS' + str(CDS_number) + '\n'+translation[0]  + '\n'                                                                       
                                FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name +'_'+'NO protein_id'+'_' + 'CDS' + str(CDS_number)+ '\n'                                                                       
                                                    
                                #FASTA_id_withEC_str ='>'+ EC_list_str +'_'+ gene_str +'_'+ product_str + ' ' + 'NO protein_id'+ ' ' + 'CDS' + str(CDS_number)+'\n'+translation[0]
                            else:
                                #FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name +'_'+'NO protein_id' +' '+ gene_str  + ' ' + EC_list_discription + ' ' + 'CDS' + str(CDS_number) + '\n'+translation[0]  + '\n'                                                                       
                                FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name +'_'+'NO protein_id'  + '\n'+translation[0]  + '\n'                                                                       
                                
                                #FASTA_id_withEC_str ='>'+ EC_list_str +'_'+ gene_str +'_'+ product_str + ' ' + 'NO protein_id'+ ' ' + 'CDS' + str(CDS_number)+'\n'+translation[0]
                            
                            
                            with open(out_fasta, 'a') as f:                           
                                print(FASTA_id_withEC_str, file = f)
                            FASTA_id_withEC.append(FASTA_id_withEC_str)
                    
                    else:
                        # To eliminate EC number(s) which includes an/other EC number(s) which is more specific
                        #if len(EC_list) >= 2:
                        
                        EC_list_nodot = [i.replace('.','') for i in EC_list]
                        EC_list_nohyphen = [i.replace('-','') for i in EC_list_nodot]
                        EC_list_nohyphen_int = [int(i) for i in EC_list_nohyphen]
                        EC_list_nohyphen_df = pd.DataFrame({'EC_num': EC_list, 'EC_num_no': EC_list_nohyphen_int})
                        
                        EC_list_nohyphen_df_sort = EC_list_nohyphen_df.sort_values('EC_num_no')               
                        EC_list_nohyphen_df_sort['EC_num_no'] = [str(i) for i in EC_list_nohyphen_df_sort['EC_num_no']]
                        EC_list_nohyphen_df_sort = EC_list_nohyphen_df_sort.reset_index(drop=True)
    
    
                        for num in EC_list_nohyphen_df_sort['EC_num_no']: 
                            for number in EC_list_nohyphen_df_sort['EC_num_no']:
                                if num == number:
                                    continue
                                elif num.startswith(number) == True:  
                                    EC_list_nohyphen_df_sort = EC_list_nohyphen_df_sort.drop(index = 0)
                                    EC_list_nohyphen_df_sort = EC_list_nohyphen_df_sort.reset_index(drop=True)
                                    continue
                                else:
                                    continue                 
                        EC_list = EC_list_nohyphen_df_sort['EC_num'] 
                        EC_list = EC_list.values.tolist()
    
                        if 'protein_id' in feature.qualifiers:
                            if 'gene' in feature.qualifiers:
                               gene = feature.qualifiers['gene'] 
                               gene_str = ''.join(gene)
                            else:
                                gene_str = 'No_gene'
                                
                            protein_id = feature.qualifiers['protein_id']
                            product = feature.qualifiers['product'] 
                            
                            EC_list = natsorted(EC_list)

                            if len(EC_list) >= 2 :
                                EC_list_3 = [ i[:5] for i in EC_list]
                                EC_list_3_unique = list(dict.fromkeys(EC_list_3))
                                c = collections.Counter(EC_list_3)
                                EC_list_3_multiple = [i[0] for i in c.items() if i[1] >= 2]
                            
                                if len(EC_list) != len(EC_list_3_unique):
                                    new_EC_list = EC_list
                                    for ECnumber_3 in EC_list_3_multiple:
                                        new_EC_list = [s for s in new_EC_list if (ECnumber_3 in s) == False]
                                        ECnumber_3_multiple = ECnumber_3 + '.X'
                                        new_EC_list.append(ECnumber_3_multiple)
                                        new_EC_list = natsorted(new_EC_list)
                                
                                    if len(new_EC_list) > 3:
                                        new_EC_list = new_EC_list[:3]
                                        EC_list_str = '+' + '_'.join(new_EC_list) 
                                    else:
                                        EC_list_str = '_'.join(new_EC_list)

                                else:
                                    if len(EC_list) < 3:
                                        EC_list_str = '_'.join(EC_list)
                                    else:
                                        new_EC_list = EC_list[:3]
                                        EC_list_str = '+' + '_'.join(new_EC_list)
                            else:
                                EC_list_str = '_'.join(EC_list)
                                                            
                            EC_list_discription = ' '.join(EC_list)   
                            protein_id_str = ' '.join(protein_id)
                            product_str = ''.join(product)                        


                            if  'translation' in feature.qualifiers:
                                translation =  feature.qualifiers['translation'] 
                            else:
                                continue
                            
                            if len(EC_list_str) < 10:  
                                #FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name + '_' + protein_id_str +'_' + 'CDS' + str(CDS_number)+ ' ' + EC_list_discription +' '+ gene_str + ' ' + 'CDS' + str(CDS_number)  + '\n'+translation[0]  + '\n'              
                                FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name + '_' + protein_id_str +'_' + 'CDS' + str(CDS_number)+ '\n'+translation[0]  + '\n'              
                                
                                #FASTA_id_withEC_str ='>'+ EC_list_str+'_'+ gene_str +'_'+ product_str + ' ' + protein_id_str + ' ' + 'CDS' + str(CDS_number)+'\n'+translation[0]        
                            else:
                                 #FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name + '_' + protein_id_str + ' ' + EC_list_discription +' '+ gene_str + ' ' + 'CDS' + str(CDS_number)  + '\n'+translation[0]  + '\n'              
                                 FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name + '_' + protein_id_str + '\n'              
                                                     
                            with open(out_fasta, 'a') as f:      
                                print(FASTA_id_withEC_str, file = f)
                            FASTA_id_withEC.append(FASTA_id_withEC_str)
    
                        else:
                            if 'gene' in feature.qualifiers:
                               gene = feature.qualifiers['gene'] 
                               gene_str = ''.join(gene)
                            else:
                                gene_str = 'No_gene'
                                
                            protein_id = [None]
                            product = feature.qualifiers['product']  
                            
                            EC_list = natsorted(EC_list)

                            if len(EC_list) >= 2 :
                                EC_list_3 = [ i[:5] for i in EC_list]
                                EC_list_3_unique = list(dict.fromkeys(EC_list_3))
                                c = collections.Counter(EC_list_3)
                                EC_list_3_multiple = [i[0] for i in c.items() if i[1] >= 2]
                            
                                if len(EC_list) != len(EC_list_3_unique):
                                    new_EC_list = EC_list
                                    for ECnumber_3 in EC_list_3_multiple:
                                        new_EC_list = [s for s in new_EC_list if (ECnumber_3 in s) == False]
                                        ECnumber_3_multiple = ECnumber_3 + '.X'
                                        new_EC_list.append(ECnumber_3_multiple)
                                        new_EC_list = natsorted(new_EC_list)
                                
                                    if len(new_EC_list) > 3:
                                        new_EC_list = new_EC_list[:3]
                                        EC_list_str = '+' + '_'.join(new_EC_list) 
                                    else:
                                        EC_list_str = '_'.join(new_EC_list)

                                else:
                                    if len(EC_list) < 3:
                                        EC_list_str = '_'.join(EC_list)
                                    else:
                                        new_EC_list = EC_list[:3]
                                        EC_list_str = '+' + '_'.join(new_EC_list)
                            else:
                                EC_list_str = '_'.join(EC_list)
                                    
                            EC_list_discription = ' '.join(EC_list)                 
                            product_str = ''.join(product)
                            
                            if  'translation' in feature.qualifiers:
                                translation =  feature.qualifiers['translation'] 
                            else:
                                continue
                            if len(EC_list_str) < 10:                                                        
                                #FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name + '_' + 'NO protein_id' +'_' + 'CDS' + str(CDS_number)+ ' ' + gene_str + ' ' + EC_list_discription + ' ' + 'CDS' + str(CDS_number) + '\n'+translation[0]+ '\n'                                                                                                    
                                FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name + '_' + 'NO protein_id' +'_' + 'CDS' + str(CDS_number) + '\n'+translation[0]+ '\n'                                                                                                    
                                
                                #FASTA_id_withEC_str ='>'+ EC_list_str +'_'+ gene_str +'_'+ product_str + ' ' + 'NO protein_id'+ ' ' + 'CDS' + str(CDS_number)+'\n'+translation[0]
                            else:
                                #FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name + '_' + 'NO protein_id' + ' ' + gene_str + ' ' + EC_list_discription + ' ' + 'CDS' + str(CDS_number) + '\n'+translation[0]+ '\n'                                                                                                    
                                FASTA_id_withEC_str ='>' + EC_list_str +'_'+ species_name + '_' + 'NO protein_id' + '\n'+translation[0]+ '\n'                                                                                                    
                                                     
                            with open(out_fasta, 'a') as f:  
                                print(FASTA_id_withEC_str, file = f)
                            FASTA_id_withEC.append(FASTA_id_withEC_str)                                     
        
                else:
                    CDS_withoutproduct_number += 1
                    
                    if 'gene' in feature.qualifiers:
                        if 'protein_id' in feature.qualifiers:
                            protein_id = feature.qualifiers['protein_id']
                            gene = feature.qualifiers['gene'] 
                            product = feature.qualifiers['product'] 
                            
                            protein_id_str = '_'.join(protein_id)
                            product_str = ''.join(product) 
                            product_str = product_str.replace(' ','_')
                            gene_str = '_'.join(gene) 
    
                            if  'translation' in feature.qualifiers:
                                translation =  feature.qualifiers['translation'] 
                            else:
                                continue
                                                        
                            #FASTA_id_withEC_str ='>NoEC_'+ species_name + '_' + gene_str + '_' + protein_id_str +'_' + 'CDS' + str(CDS_number)+' ' + product_str + '\n'+translation[0]  + '\n'       
                            FASTA_id_withEC_str ='>NoEC_'+ species_name + '_' + gene_str + '_' + protein_id_str +'_' + 'CDS' + str(CDS_number)+ '\n'+translation[0]  + '\n'       
                            
                            #FASTA_id_withEC_str ='>'+ EC_list_str+'_'+ gene_str +'_'+ product_str + ' ' + protein_id_str + ' ' + 'CDS' + str(CDS_number)+'\n'+translation[0]        
                            
                            with open(out_fasta, 'a') as f:      
                                print(FASTA_id_withEC_str, file = f)
                            FASTA_id_withEC.append(FASTA_id_withEC_str)
    
                        else:
                            protein_id = [None]
                            gene = feature.qualifiers['gene'] 
                            product = feature.qualifiers['product']  
                                            
                            product_str = ''.join(product)
                            gene_str = ''.join(gene)
                            
                            if  'translation' in feature.qualifiers:
                                translation =  feature.qualifiers['translation'] 
                            else:
                                continue
                            
                            FASTA_id_withEC_str ='>NoEC_'+ species_name + '_' + gene_str + '_' + '_NO protein_id'+'_' + 'CDS' + str(CDS_number) + '\n'+translation[0]  + '\n'                                                          
                           #FASTA_id_withEC_str ='>'+ EC_list_str +'_'+ gene_str +'_'+ product_str + ' ' + 'NO protein_id'+ ' ' + 'CDS' + str(CDS_number)+'\n'+translation[0]
                           
                            with open(out_fasta, 'a') as f:  
                                print(FASTA_id_withEC_str, file = f)
                            FASTA_id_withEC.append(FASTA_id_withEC_str)
                    else:
                        if 'protein_id' in feature.qualifiers:
                            protein_id = feature.qualifiers['protein_id']
                            product = feature.qualifiers['product'] 
                            
                            protein_id_str = '_'.join(protein_id)
                            product_str = ''.join(product) 
                            product_str = product_str.replace(' ','_')
                            gene_str = 'No_gene' 
    
                            if  'translation' in feature.qualifiers:
                                translation =  feature.qualifiers['translation'] 
                            else:
                                continue
                                                        
                            #FASTA_id_withEC_str ='>NoEC_'+ species_name + '_' + gene_str + '_' + protein_id_str +'_' + 'CDS' + str(CDS_number) +' ' + product_str+ '\n'+translation[0]  + '\n'       
                            FASTA_id_withEC_str ='>NoEC_'+ species_name + '_' + gene_str + '_' + protein_id_str +'_' + 'CDS' + str(CDS_number) + '\n'+translation[0]  + '\n'       
                            
                            #FASTA_id_withEC_str ='>'+ EC_list_str+'_'+ gene_str +'_'+ product_str + ' ' + protein_id_str + ' ' + 'CDS' + str(CDS_number)+'\n'+translation[0]        
                            
                            with open(out_fasta, 'a') as f:      
                                print(FASTA_id_withEC_str, file = f)
                            FASTA_id_withEC.append(FASTA_id_withEC_str)
    
                        else:
                            protein_id = [None]
                            product = feature.qualifiers['product']  
                                            
                            product_str = ''.join(product)
                            gene_str = 'No_gene'
                            
                            if  'translation' in feature.qualifiers:
                                translation =  feature.qualifiers['translation'] 
                            else:
                                continue
                            
                            #FASTA_id_withEC_str ='>NoEC_'+ species_name + '_' + gene_str + '_' + '_NO protein_id'  +'_' + 'CDS' + str(CDS_number)+' ' + product_str+ '\n'+translation[0]  + '\n'                                                          
                            FASTA_id_withEC_str ='>NoEC_'+ species_name + '_' + gene_str + '_' + '_NO protein_id'  +'_' + 'CDS' + str(CDS_number)+ '\n'+translation[0]  + '\n'                                                          
                          
                            #FASTA_id_withEC_str ='>'+ EC_list_str +'_'+ gene_str +'_'+ product_str + ' ' + 'NO protein_id'+ ' ' + 'CDS' + str(CDS_number)+'\n'+translation[0]
                           
                            with open(out_fasta, 'a') as f:  
                                print(FASTA_id_withEC_str, file = f)
                            FASTA_id_withEC.append(FASTA_id_withEC_str)                            
                    
                        continue    




if __name__ == '__main__':

    path = 'C:/work_master/全ゲノムデータ_blast2回目'   
    files = os.listdir(path)
    
    for file in files:
        gb_file = file
        CDS_count = 0
        #gb_file = 'multiECnumber.gb'
        for record in SeqIO.parse(gb_file, "genbank"):       
            for feature in record.features:
                if feature.type == "CDS":
                    CDS_count += 1
                else:
                    continue
    
        # output file name construction
       
        if CDS_count == 1:
            out_fasta = gb_file + '_all-CDS_ECnum_1.fa'
        elif CDS_count >=2:
            out_fasta = gb_file + '_all-CDS_ECnum_' + str(CDS_count) + '.mfa'
        else: # por si acaso...
            out_fasta = gb_file + '_NO-CDS-FOUND.txt'
            check = 0
    
    
        make_fasta(gb_file,out_fasta,CDS_count)

    
