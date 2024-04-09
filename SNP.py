"""
Created on Thu Jun 23 21:00:00 2022

@author: Hugo Leblond - Master 2 AVR / Telecom Nancy / Université de Lorraine
@contact: hugo.leblond@telecomnancy.eu - hugo.leblond7@gmail.com
"""
################################## Main goal ##################################
# The purpose of these functions is to identify the parental origin of chromosomal
# regions in a recombinant. The algorithm first assigns each SNP to one of the 
# two parents (P1 or P2). It then defines regions assigned to P1 or P2 by 
# applying a threshold of minimum number of SNPs (typically 2), then calculates 
# their size relative to the parental genomes.
###############################################################################


import pandas as pd
import xlsxwriter
from datetime import datetime


def SNP_parent(nom_fichier_excel):
    """
    Assign each SNP to a parent with its position in P1, P2 and R.
    
    Parameters
    ----------
    nom_fichier_excel : String,
        SNP file name from MAUVE converted to Excel format.

    Returns
    -------
    liste_appartenance_SNP : list,
        List containing the labels of the SNPs (P1 or P2) in string format.
    liste_position_SNP_P1 : list,
        List of position of each SNP related to the P1 genome.
    liste_position_SNP_P2 : list,
        List of position of each SNP related to the P2 genome.
    liste_position_SNP_R : list,
        List of position of each SNP related to the Recombinant consensus genome (see material and methods in the related article).

    """
    mon_fichier = pd.read_excel(nom_fichier_excel) 
    SNP_pattern = list(mon_fichier['SNP pattern'])  
    
    positions_SNP_P1 = list(mon_fichier['sequence_1_PosInContg'])
    positions_SNP_P2 = list(mon_fichier['sequence_2_PosInContg'])
    positions_SNP_R = list(mon_fichier['sequence_3_PosInContg'])
    
    liste_appartenance_SNP = []                    
    liste_position_SNP_P1 = []
    liste_position_SNP_P2 = []
    liste_position_SNP_R = []
    
    compt = 0
    liste_SNP_pattern = [] 
               
    for snp_pattern in SNP_pattern :                
        liste_SNP_pattern.append(str(snp_pattern))  
        
    for snp in liste_SNP_pattern :                  
        
        if snp[2].lower() == 'n' or snp[1].lower() == 'n' or snp[0].lower() == 'n' :  
            compt+=1
        
        elif snp[2].lower() == snp[0].lower():         
            if positions_SNP_P1[compt] != 0 and positions_SNP_P2[compt] != 0 :
                liste_appartenance_SNP.append("P1")     
                liste_position_SNP_P1.append(positions_SNP_P1[compt])
                liste_position_SNP_P2.append(positions_SNP_P2[compt])
                liste_position_SNP_R.append(positions_SNP_R[compt])
                
            compt+=1
            
        elif snp[2].lower() == snp[1].lower() :
            if positions_SNP_P1[compt] != 0 and positions_SNP_P2[compt] != 0 :
                liste_appartenance_SNP.append("P2")
                liste_position_SNP_P1.append(positions_SNP_P1[compt])
                liste_position_SNP_P2.append(positions_SNP_P2[compt])
                liste_position_SNP_R.append(positions_SNP_R[compt])
            compt+=1
        
        else :                                      
            compt+=1
        
    return liste_appartenance_SNP, liste_position_SNP_P1, liste_position_SNP_P2, liste_position_SNP_R   



def SNP_chaine(liste_appartenance_SNP, base=1, nombre_suite_P1=2, nombre_suite_P2=2):
    """
    Define regions of identical successive SNP labels.

    Parameters
    ----------
    liste_appartenance_SNP : list,
        List containing the labels of the SNPs (P1 or P2) in string format.
    base : int, by default the receptor genome (1 or 2 according to P1 and P2 parents), automatically selected in the main file.
    nombre_suite_P1 : int, optional
        Minimum number of successive SNP belonging to the same parent P1. The default is 2.
    nombre_suite_P2 : int, optional
        Minimum number of successive SNP belonging to the same parent P2. The default is 2.

    Returns
    -------
    chaine_SNP : list,
        List containing the labels of the SNPs (P1 or P2) in string format. Filter positions under the threshold defined by nombre_suite_P1 and nombre_suite_P2.

    """
    chaine_SNP = []
    init_boucle = 0
    compt_suite = 0
    switch1 = "P1"
    switch2 = "P2"
    switch = ""
    base = "P"+str(base)
    
    if base == "P1":
        switch = switch2
    else :
        switch = switch1
        
    while init_boucle <= len(liste_appartenance_SNP)-1 :   
        if switch == switch1 :                              
            if liste_appartenance_SNP[init_boucle] == switch2 : 
                compt_suite += 1                                
                chaine_SNP.append(1)                            
                init_boucle+=1
            elif liste_appartenance_SNP[init_boucle]==switch1 : 
                compt_suite = 0
                chaine_SNP.append(1)
                init_boucle+=1
            
            
            if compt_suite == nombre_suite_P2 :                    
                for i in range(nombre_suite_P2):                   
                    chaine_SNP[init_boucle - i - 1]=2
                compt_suite = 0                                
                switch = switch2                              
                
        else :
            if liste_appartenance_SNP[init_boucle] == switch1 :
                compt_suite += 1                                
                chaine_SNP.append(2)                            
                init_boucle+=1
            elif liste_appartenance_SNP[init_boucle]==switch2 : 
                compt_suite = 0
                chaine_SNP.append(2)                            
                init_boucle+=1
               
            
            if compt_suite == nombre_suite_P1 :                    
                for i in range(nombre_suite_P1):                   
                    chaine_SNP[init_boucle - i - 1]=1
                compt_suite = 0                                
                switch = switch1                               
            
    return chaine_SNP

def compteur_snp(chaine_snp):
    """
    Counts the number of SNPs in each segment.

    Parameters:
    - chaine_snp (list): List representing the SNP sequence.

    Returns:
    - liste_taille_snp (list): List containing the number of SNPs in each segment.
    """
    compt_p1 = 0
    compt_p2 = 0
    liste_taille_snp = []

    for i in range(len(chaine_snp)-1):

        if i+1 == len(chaine_snp)-1:
            if compt_p1 != 0:
                compt_p1 += 1
                liste_taille_snp.append(compt_p1)
            else:
                compt_p2 += 1
                liste_taille_snp.append(compt_p2)

        elif chaine_snp[i] == 2 and chaine_snp[i+1] == 1:
            compt_p2 += 1
            liste_taille_snp.append(compt_p2)
            compt_p2 = 0
        elif chaine_snp[i] == 1 and chaine_snp[i+1] == 2:
            compt_p1+=1
            liste_taille_snp.append(compt_p1)
            compt_p1 = 0

        elif chaine_snp[i] == 1:
            compt_p1 += 1
        elif chaine_snp[i] == 2:
            compt_p2 += 1
    return liste_taille_snp



#Convertit les chaines du type [P1,P1,P1 ... P2,P2,P2... P1,P1,P1] avec leur position en segments (ex : P1, debutP1,finP1)
def chaine_borne(chaine_SNP, liste_P1, liste_P2, liste_R, donneur=2):
    """
    Take a SNP list and lists of position of each SNP related to the P1, P2 and R genome. 
    Return lists of positions (start/end) of stretches of SNP related to P1, P2 and R.

    Parameters
    ----------
    chaine_SNP : list,
        List containing the labels of the SNPs (P1 or P2) in string format.
    liste_P1 : list,
        List of position of each SNP related to the P1 genome.
    liste_P2 : list,
        List of position of each SNP related to the P2 genome.
    liste_R : list,
        List of position of each SNP related to the R genome.
    donneur : int, optional
        by default the donor genome (1 or 2 according to P1 and P2 parents), automatically selected in the main file.

    Returns
    -------
    chaine_borne_SNP : list,
        List containing the labels of the SNPs (P1 or P2) in string format.
    chaine_borne_P1_start : list,
        List of start positions of SNP stretches relative to P1.
    chaine_borne_P2_start : list,
        List of start positions of SNP stretches relative to P2.
    chaine_borne_R_start : list,
        List of start positions of SNP stretches relative to R.
    chaine_borne_P1_end : list,
        List of end positions of SNP stretches relative to P1.
    chaine_borne_P2_end : list,
        List of end positions of SNP stretches relative to P2.
    chaine_borne_R_end : list,
        List of end positions of SNP stretches relative to R.

    """
    chaine_borne_SNP = []
    chaine_borne_P1_start = []
    chaine_borne_P2_start = []
    chaine_borne_R_start = []
    
    chaine_borne_P1_end = []
    chaine_borne_P2_end = []
    chaine_borne_R_end = []
    
    init_boucle = 0  
    
    while init_boucle <= len(chaine_SNP)-1 :
        
        if init_boucle == 0 :
            chaine_borne_SNP.append(chaine_SNP[init_boucle])    
            chaine_borne_P1_start.append(liste_P1[init_boucle])
            chaine_borne_P2_start.append(liste_P2[init_boucle])
            chaine_borne_R_start.append(liste_R[init_boucle])
            init_boucle+=1
        
     
        elif init_boucle == len(chaine_SNP)-1 :
            chaine_borne_P1_end.append(liste_P1[init_boucle])
            chaine_borne_P2_end.append(liste_P2[init_boucle])
            chaine_borne_R_end.append(liste_R[init_boucle])
            init_boucle+=1
            
 
        elif chaine_SNP[init_boucle]!=chaine_SNP[init_boucle+1] and init_boucle+1 != len(chaine_SNP)-1:
            if chaine_SNP[init_boucle] == donneur :
                
                chaine_borne_P1_end.append(liste_P1[init_boucle]) 
                chaine_borne_P2_end.append(liste_P2[init_boucle])
                chaine_borne_R_end.append(liste_R[init_boucle])
                
                chaine_borne_SNP.append(chaine_SNP[init_boucle+1])
                chaine_borne_P1_start.append(liste_P1[init_boucle]+1)
                chaine_borne_P2_start.append(liste_P2[init_boucle]+1)
                chaine_borne_R_start.append(liste_R[init_boucle]+1)
            else :
                chaine_borne_P1_end.append(liste_P1[init_boucle+1]-1)
                chaine_borne_P2_end.append(liste_P2[init_boucle+1]-1)
                chaine_borne_R_end.append(liste_R[init_boucle+1]-1)
                
                chaine_borne_SNP.append(chaine_SNP[init_boucle+1])
                chaine_borne_P1_start.append(liste_P1[init_boucle+1])
                chaine_borne_P2_start.append(liste_P2[init_boucle+1])
                chaine_borne_R_start.append(liste_R[init_boucle+1])
            
            init_boucle+=1
        
        else :
            init_boucle+=1
        
        
                                   
    return chaine_borne_SNP, chaine_borne_P1_start, chaine_borne_P2_start, chaine_borne_R_start, chaine_borne_P1_end, chaine_borne_P2_end, chaine_borne_R_end

def calcul_des_sommes_rapport_recombinant(chaine_borne_SNP, chaine_borne_Receptor_start, chaine_borne_Receptor_end, borne_sup =1000):
    """

    Parameters
    ----------
    chaine_borne_SNP : list,
        List containing the labels of the SNPs (P1 or P2) in string format.
    chaine_borne_Receptor_start : list,
        DESCRIPTION.
    chaine_borne_Receptor_end : list,
        List of start positions of SNP stretches relative to the receptor (P1 or P2). Automatically selected in the main file.
    borne_sup : int, optional
        Compute the sum of regions of length exceeding the threshold. The default is 1000. Not displayed in the excel output.

    Returns
    -------
    somme_P1_liste : list,
        List containing the total length of P1 regions.
    somme_P2_liste : list,
        List containing the total length of P2 regions.
    min_P1_P2_liste : String,
        Label of the parent with the minimum total length.
    somme_sup_liste : list,
        DESCRIPTION.
    compt : list,
        List containing the number of fragments over the threshold length borne_sup. Not displayed in the output excel file.
    compt_P1 : list,
        List containing the number of fragments of label P1.
    compt_P2 : list,
        List containing the number of fragments of label P2.

    """
    somme_P1 = 0
    somme_P2 = 0
    somme_sup = 0
    
    compt_P1 = [0]
    compt_P2 = [0]
    compt = [0]
    
    somme_P1_liste = [] 
    somme_P2_liste = []
    min_P1_P2_liste = []
    somme_sup_liste = []
    
    for i in range(len(chaine_borne_SNP)):
        if chaine_borne_SNP[i] == 1 :
            compt_P1[0] +=1
            somme_P1 += chaine_borne_Receptor_end[i] - chaine_borne_Receptor_start[i]
            
        elif chaine_borne_SNP[i] == 2 : 
            compt_P2[0] +=1
            somme_P2 += chaine_borne_Receptor_end[i] - chaine_borne_Receptor_start[i]
                
    somme_P1_liste.append(somme_P1)
    somme_P2_liste.append(somme_P2)
    
    # Partie pour le calcul des sommes bornées
    if min(somme_P1, somme_P2) == somme_P1 : 
        min_P1_P2_liste.append("P1")
        for i in range(len(chaine_borne_SNP)):
            if chaine_borne_SNP[i] == 1 and chaine_borne_Receptor_end[i] - chaine_borne_Receptor_start[i] > borne_sup :
                compt[0]+=1
                somme_sup += chaine_borne_Receptor_end[i] - chaine_borne_Receptor_start[i]
    elif min(somme_P1, somme_P2) == somme_P2 : 
        min_P1_P2_liste.append("P2")
        for i in range(len(chaine_borne_SNP)):
            if chaine_borne_SNP[i] == 2 and chaine_borne_Receptor_end[i] - chaine_borne_Receptor_start[i] > borne_sup :
                somme_sup += chaine_borne_Receptor_end[i] - chaine_borne_Receptor_start[i]
                compt[0]+=1

    
    somme_sup_liste.append(somme_sup)
    
    return somme_P1_liste, somme_P2_liste, min_P1_P2_liste, somme_sup_liste, compt, compt_P1, compt_P2

def calcul_des_sommes_rapport_donneur(chaine_borne_SNP, chaine_borne_donneur_start, chaine_borne_donneur_end):
    """

    Parameters
    ----------
    chaine_borne_SNP :  list,
        List containing the labels of the SNPs (P1 or P2) in string format.
    chaine_borne_donneur_start : list,
        List of start positions of SNP stretches relative to the donor (P1 or P2). Automatically selected in the main file.
    chaine_borne_donneur_end : list,
        List of end positions of SNP stretches relative to the donor (P1 or P2). Automatically selected in the main file.

    Returns
    -------
    somme_P1_liste : list,
        List containing the total length of P1 regions.
    somme_P2_liste : list,
        List containing the total length of P2 regions.

    """
    somme_P1 = 0
    somme_P2 = 0
    
    somme_P1_liste = []
    somme_P2_liste = []
    
    for i in range(len(chaine_borne_SNP)):
        if chaine_borne_SNP[i] == 1 :
            somme_P1 += chaine_borne_donneur_end[i] - chaine_borne_donneur_start[i]
        elif chaine_borne_SNP[i] == 2 : 
            somme_P2 += chaine_borne_donneur_end[i] - chaine_borne_donneur_start[i]
                
    somme_P1_liste.append(somme_P1) 
    somme_P2_liste.append(somme_P2)
    return somme_P1_liste, somme_P2_liste


def creer_excel(nom_fichier_entree, chaine_borne_SNP, liste_debut, liste_fin, liste_P1_debut, liste_P1_fin, liste_P2_debut, liste_P2_fin,chaine_borne_cote_donneur,liste_cote_donneur_start, liste_cote_donneur_end, somme_P1, somme_P2, min_P1_P2, compt_P1, compt_P2, somme_donneur_indices_donneur, liste_somme):
    """
    

    Parameters
    ----------
    nom_fichier_entree : String,
        SNP file name from MAUVE converted to Excel format.
    chaine_borne_SNP : list,
        List containing the labels of the SNPs (P1 or P2) in string format.
    liste_debut : list,
        List of start positions of SNP stretches relative to R.
    liste_fin : list,
        List of end positions of SNP stretches relative to R.
    liste_P1_debut : list,
        List of start positions of SNP stretches relative to P1.
    liste_P1_fin : list,
        List of end positions of SNP stretches relative to P1.
    liste_P2_debut : list,
        List of start positions of SNP stretches relative to P2.
    liste_P2_fin : list,
        List of end positions of SNP stretches relative to P2.
    chaine_borne_cote_donneur : list,
        List containing the labels of the SNPs (P1 or P2) in string format sorted relative to the donor.
    liste_cote_donneur_start : list,
        List containing the start positions of the SNP regions relative to the donor.
    liste_cote_donneur_end : list,
        List containing the end positions of the SNP regions relative to the donor.
    somme_P1 : list,
        List containing the total length of P1 regions.
    somme_P2 : list,
        List containing the total length of P2 regions.
    min_P1_P2 : String,
        Label of the parent with the minimum total length.
    compt_P1 : list,
        List containing the number of fragments of label P1.
    compt_P2 : list,
        List containing the number of fragments of label P2.
    somme_donneur_indices_donneur : list,
        Total length of donor fragments relative to the donor sorted list. 

    Returns
    -------
    Excel file containing all the results.

    """
    now = datetime.now()
    nom_fichier = nom_fichier_entree[: (len(nom_fichier_entree) - 5)] +"_sortie"  + now.strftime("_%d_%m_%Y %H_%M") + ".xlsx"
    xlsxwriter.Workbook(nom_fichier)
    
    donnee = pd.DataFrame({"Pattern" : chaine_borne_SNP,
                  "Start" : liste_debut, 
                  "End" : liste_fin})
    
    donnee_P1_P2_cote_recepteur = pd.DataFrame({"Pattern" : chaine_borne_SNP,
                                      "Nombre de SNP" : liste_somme,
                                      "Start P1" : liste_P1_debut,
                                      "End P1" : liste_P1_fin,
                                      "Start P2" : liste_P2_debut,
                                      "End P2" : liste_P2_fin,
                                      })
    sommes_cote_recepteur = pd.DataFrame({
                                      "somme P1" : somme_P1,
                                      "somme P2" : somme_P2, 
                                      "min P1 P2" : min_P1_P2,
                                      "nombre de fragments P1" : compt_P1,
                                      "nombre de fragments P2" : compt_P2,})
    
    donnee_P2_cote_donneur = pd.DataFrame({"Pattern" : chaine_borne_cote_donneur,
                  "Start P2" : liste_cote_donneur_start,
                  "End P2" : liste_cote_donneur_end,})
    somme_cote_donneur = pd.DataFrame({
                  "Somme du donneur" : somme_donneur_indices_donneur,})

    with pd.ExcelWriter(nom_fichier) as writer:
        donnee.to_excel(writer, sheet_name='donnees R')
        donnee_P1_P2_cote_recepteur.to_excel(writer, sheet_name="Cotes recepteur")
        sommes_cote_recepteur.to_excel(writer, sheet_name="Somme recepteur")
        donnee_P2_cote_donneur.to_excel(writer, sheet_name='Cotes donneur')
        somme_cote_donneur.to_excel(writer, sheet_name='Somme donneur')

    
    return 0




