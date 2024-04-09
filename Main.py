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

from SNP import SNP_parent, SNP_chaine, chaine_borne, calcul_des_sommes_rapport_recombinant, calcul_des_sommes_rapport_donneur, creer_excel, compteur_snp


# Insert the input files names with your exact excel files names (SNP MAUVE output files).
# File sorted by positions on either P1
Nom_fichier_tri_P1 = "SNP_R11_consRLB1_9_triéP1.xlsm"
# and P2 parents.
Nom_fichier_tri_P2 = "SNP_R11_consRLB1_9_triéP2.xlsm"


def utiliser_fonctions(Nom_fichier_tri_P1, Nom_fichier_tri_P2, nb_suite_P1, nb_suite_P2):
    """
    Automatically select the donor and receptor parent/genome and uses the SNP package functions to create the output excel file.

    Parameters
    ----------
    Nom_fichier_tri_P1 : String, 
        Name of the excel file sorted by positions on P1 genome.
    Nom_fichier_tri_P2 : String, 
        Name of the excel file sorted by positions on P2 genome.
    nb_suite_P1 : int,
        Minimum number of successive SNP belonging to the same parent P1. The default is 2.
    nb_suite_P2 : int,
        Minimum number of successive SNP belonging to the same parent P2. The default is 2.

    Returns
    -------
    Excel file containing all the results.

    """
    Recepteur = 1   # Le prédicat de base est que P1 est le recepteur et P2 le donneur. (Ces variables sont modifiées automatiquement)
    Donneur = 2
    #Par rapport à P1
    liste_parent, liste_P1, liste_p2, liste_R = SNP_parent(Nom_fichier_tri_P1)
    chaine_SNP = SNP_chaine(liste_parent, Recepteur, nb_suite_P1, nb_suite_P2) # variables : 1 --> base (ex: P1), 3 --> nb de snp d'affilé pour passer de P1 a P2 et réciproquement
    liste_taille_snp = compteur_snp(chaine_SNP)
    chaine_borne_SNP, chaine_borne_P1_start, chaine_borne_P2_start, chaine_borne_R_start, chaine_borne_P1_end, chaine_borne_P2_end, chaine_borne_R_end = chaine_borne(chaine_SNP, liste_P1, liste_p2, liste_R, Donneur) # Variable 2 --> donneur (valeur possible 1 ou 2)
    somme_P1, somme_P2, min_P1_P2, somme_sup, compt, compt_P1, compt_P2 = calcul_des_sommes_rapport_recombinant(chaine_borne_SNP, chaine_borne_P1_start, chaine_borne_P1_end, 1000)
    
    #Calcul de somme sur P2 trié
    liste_parent_, liste_P1_, liste_p2_, liste_R_ = SNP_parent(Nom_fichier_tri_P2)
    chaine_SNP_ = SNP_chaine(liste_parent_, Recepteur, nb_suite_P1, nb_suite_P2)
    # liste_taille_snp = compteur_snp(chaine_SNP_)
    chaine_borne_SNP_, _, chaine_borne_P2_start_, _, _, chaine_borne_P2_end_, _ = chaine_borne(chaine_SNP_, liste_P1_, liste_p2_, liste_R_, Donneur)
    somme_P1_donneur, somme_P2_donneur = calcul_des_sommes_rapport_donneur(chaine_borne_SNP_, chaine_borne_P2_start_, chaine_borne_P2_end_)
    if min_P1_P2[0] == "P2":
        creer_excel(Nom_fichier_tri_P1, chaine_borne_SNP, chaine_borne_R_start, chaine_borne_R_end, chaine_borne_P1_start, chaine_borne_P1_end, chaine_borne_P2_start, chaine_borne_P2_end,chaine_borne_SNP_,chaine_borne_P2_start_,chaine_borne_P2_end_, somme_P1, somme_P2, min_P1_P2, compt_P1, compt_P2, somme_P2_donneur, liste_taille_snp)
        print("Opération terminée : P1 est le recepteur et P2 le donneur")
     
        Recepteur = 2
        Donneur = 1
        #Par rapport a P2
        liste_parent, liste_P1, liste_p2, liste_R = SNP_parent(Nom_fichier_tri_P2)
        chaine_SNP = SNP_chaine(liste_parent, Recepteur, nb_suite_P1, nb_suite_P2) # variables : 1 --> base (ex: P1), 3 --> nb de snp d'affilé pour passer de P1 a P2 et réciproquement
        chaine_borne_SNP, chaine_borne_P1_start, chaine_borne_P2_start, chaine_borne_R_start, chaine_borne_P1_end, chaine_borne_P2_end, chaine_borne_R_end = chaine_borne(chaine_SNP, liste_P1, liste_p2, liste_R, Donneur) # Variable 2 --> donneur (valeur possible 1 ou 2)
        somme_P1, somme_P2, min_P1_P2, somme_sup, compt, compt_P1, compt_P2 = calcul_des_sommes_rapport_recombinant(chaine_borne_SNP, chaine_borne_P2_start, chaine_borne_P2_end, 1000)
        
        #Calcul de somme sur P1 trié
        liste_parent_, liste_P1_, liste_p2_, liste_R_ = SNP_parent(Nom_fichier_tri_P1)
        chaine_SNP_ = SNP_chaine(liste_parent_, Recepteur, nb_suite_P1, nb_suite_P2) # variables : 1 --> base (ex: P1), 3 --> nb de snp d'affilé pour passer de P1 a P2 et réciproquement
        liste_taille_snp = compteur_snp(chaine_SNP)
        chaine_borne_SNP_, chaine_borne_P1_start_, _, _, chaine_borne_P1_end_, _, _ = chaine_borne(chaine_SNP_, liste_P1_, liste_p2_, liste_R_, Donneur)
        somme_P1_donneur, somme_P2_donneur = calcul_des_sommes_rapport_donneur(chaine_borne_SNP_, chaine_borne_P1_start_, chaine_borne_P1_end_)
        #somme_P1_donneur, _, _, _,_, _, _ = calcul_des_sommes_rapport_recombinant(chaine_borne_SNP, chaine_borne_P1_start_, chaine_borne_P1_end_, 1000)
        
        creer_excel(Nom_fichier_tri_P2, chaine_borne_SNP, chaine_borne_R_start, chaine_borne_R_end, chaine_borne_P1_start, chaine_borne_P1_end, chaine_borne_P2_start, chaine_borne_P2_end,chaine_borne_SNP_ ,chaine_borne_P1_start_ ,chaine_borne_P1_end_, somme_P1, somme_P2, min_P1_P2, compt_P1, compt_P2, somme_P1_donneur, liste_taille_snp)

        return 0
    else :
        print("Opération terminée : P2 est le recepteur et P1 le donneur")
        creer_excel(Nom_fichier_tri_P1, chaine_borne_SNP, chaine_borne_R_start, chaine_borne_R_end, chaine_borne_P1_start, chaine_borne_P1_end, chaine_borne_P2_start, chaine_borne_P2_end,chaine_borne_SNP_,chaine_borne_P2_start_,chaine_borne_P2_end_, somme_P1, somme_P2, min_P1_P2, compt_P1, compt_P2, somme_P2_donneur, liste_taille_snp)
        
        Recepteur = 2
        Donneur = 1
        #Par rapport a P2
        liste_parent, liste_P1, liste_p2, liste_R = SNP_parent(Nom_fichier_tri_P2)
        chaine_SNP = SNP_chaine(liste_parent, Recepteur, nb_suite_P1, nb_suite_P2) # variables : 1 --> base (ex: P1), 3 --> nb de snp d'affilé pour passer de P1 a P2 et réciproquement
        chaine_borne_SNP, chaine_borne_P1_start, chaine_borne_P2_start, chaine_borne_R_start, chaine_borne_P1_end, chaine_borne_P2_end, chaine_borne_R_end = chaine_borne(chaine_SNP, liste_P1, liste_p2, liste_R, Donneur) # Variable 2 --> donneur (valeur possible 1 ou 2)
        somme_P1, somme_P2, min_P1_P2, somme_sup, compt, compt_P1, compt_P2 = calcul_des_sommes_rapport_recombinant(chaine_borne_SNP, chaine_borne_P2_start, chaine_borne_P2_end, 1000)
        
        #Calcul de somme sur P1 trié
        liste_parent_, liste_P1_, liste_p2_, liste_R_ = SNP_parent(Nom_fichier_tri_P1)
        chaine_SNP_ = SNP_chaine(liste_parent_, Recepteur, nb_suite_P1, nb_suite_P2) # variables : 1 --> base (ex: P1), 3 --> nb de snp d'affilé pour passer de P1 a P2 et réciproquement
        liste_taille_snp = compteur_snp(chaine_SNP)
        chaine_borne_SNP_, chaine_borne_P1_start_, _, _, chaine_borne_P1_end_, _, _ = chaine_borne(chaine_SNP_, liste_P1_, liste_p2_, liste_R_, Donneur)
        somme_P1_donneur, somme_P2_donneur = calcul_des_sommes_rapport_donneur(chaine_borne_SNP_, chaine_borne_P1_start_, chaine_borne_P1_end_)
        #somme_P1_donneur, _, _, _,_, _, _ = calcul_des_sommes_rapport_recombinant(chaine_borne_SNP, chaine_borne_P1_start_, chaine_borne_P1_end_, 1000)
        
        creer_excel(Nom_fichier_tri_P2, chaine_borne_SNP, chaine_borne_R_start, chaine_borne_R_end, chaine_borne_P1_start, chaine_borne_P1_end, chaine_borne_P2_start, chaine_borne_P2_end,chaine_borne_SNP_ ,chaine_borne_P1_start_ ,chaine_borne_P1_end_, somme_P1, somme_P2, min_P1_P2, compt_P1, compt_P2, somme_P1_donneur, liste_taille_snp)
        
        return 0

utiliser_fonctions(Nom_fichier_tri_P1, Nom_fichier_tri_P2, 2, 2)



