import sys
from comp_sim import *
from concordance_index import CI
from sequence_alignment import *
from subprocess import check_call
from read import *
import os





#################################################################################
#process of cross ranking based methods
#train a model from bioassay1, and then predict bioassay2 with bioassay1's model
#bioAssay1= path to bioassay 1
#bioAssay2= path to bioassay 2
#model_file= path to save bioassay 1's model
#pred_file= path to save bioassay 2' prediction with bioassay 1's model
#c= parameter '-c' for SVM light to train bioassay 1's model
#t= parameter '-t' for SVM light to train bioassay 1's model
#################################################################################
def cross_ranking_method(bioAssay1, bioAssay2, model_file, pred_file, c, t):


	#check if model file and prediction files are set
        if model_file == None or pred_file == None:
                sys.exit("Under \"cross ranking\" mode, Please specify a model file path and a pred file path. Please see documentation for details.")

        #train model
        #if a model alread exists, training process will be skipped
        #Make sure yo do not
        if not os.path.isfile(model_file):
                #path to save log file
                model_log = model_file + ".log"
                #check if parameters are set properly
                #Please see documentation for explanation
                #t = type of kernel. 
                #c = constrain parameter
                #-z p: preference ranking problem. In this project, this parameter is always fixed.
                if t == None or c == None:
                        sys.exit("To train a model, Please specify -t and -c parameters. Please check documentation for details.")
                #generate a command to train a model
                cmd_train_model = "./svm_learn -z p -t " + str(t) + " -c " + str(c) + " " + bioAssay1 + " " + model_file + " > " + model_log

                #train a file
                check_call(cmd_train_model, shell=True)

	
        #predict bioAssay2 with bioAssay1's model
        #path to save prediction file
        pred_log = pred_file + ".log"
        cmd_pred = "./svm_classify " + bioAssay2 + " " + model_file + " " + pred_file  + " > " + pred_log
        check_call(cmd_pred, shell=True)






#################################################################################
#calculate similarity by average pairwise compound similarity
#################################################################################
#bioAssay1= path to bioassay 1
#bioAssay2= path to bioassay 2
def pairwise_sim(bioAssay1, bioAssay2):
	#read the bioassays first
	assay1 = read_assay_feature(bioAssay1)
        assay2 = read_assay_feature(bioAssay2)

	#initialze the sum of similarities
        totalSim = 0.0

        #calculate the similairies of all pairs of compounds
        #sum up the pairwise similarities
        for comp1 in assay1:
                for comp2 in assay2:
                        totalSim += tanimoto(comp1, comp2)

        #normalize by the number of pairs
        similarity = totalSim / (float(len(assay1)) * float(len(assay2)))

	return similarity



#################################################################################
#calculate similarity by cross ranking CI
#use bioAssay1's model to predict bioAssay2's training data sets
#calculate the CI of such training data (as ground truth) and the prediction
#################################################################################
#bioAssay1= path to bioassay 1
#bioAssay2= path to bioassay 2
#model_file= path to save bioassay 1's model
#pred_file= path to save bioassay 2' prediction with bioassay 1's model
#c= parameter '-c' for SVM light to train bioassay 1's model
#t= parameter '-t' for SVM light to train bioassay 1's model
def xci_sim(bioAssay1, bioAssay2, model_file, pred_file, c, t):
	#call cross ranking based method to train bioAssay1's model
	#and predict bioAssay2 with bioAssay1's model
	cross_ranking_method(bioAssay1, bioAssay2, model_file, pred_file, c, t)	

        #calcualte concordance index for bioAssay2 and its prediction with bioAssay1's model
        similarity = CI(bioAssay2, pred_file)

	return similarity


#calculate similarity by profiling based sequence alignment by compound ID
#use bioAssay1 and bioassay2 to align
#bioAssay1= path to bioAssay1
#bioAssay2= path to bioAssay2
#match_score= user-defined score for a matched step
#mismatch_score= user-defined score for a mismatched step
#gap_score= user-defined score for a gap step
#disc = (bool) if True, alignment score is discounted by position
def p_cid_sim(bioAssay1, bioAssay2, match_score, mismatch_score, gap_score, disc):
	#read the bioassays first
	assay1 = read_assay_feature(bioAssay1)
	assay2 = read_assay_feature(bioAssay2)

	#initialize parameters
	align_mode = "cid"

	#calculate similarity
	similarity = sequence_alignment_similarity(assay1, assay2, align_mode, match_score, mismatch_score, gap_score, disc)

	return similarity
	



#calculate similarity by profiling based sequence alignment by compound similarity
#use bioAssay1 and bioassay2 to align
#bioAssay1= path to bioAssay1
#bioAssay2= path to bioAssay2
#disc = (bool) if True, alignment score is discounted by position
def p_csim_sim(bioAssay1, bioAssay2, gap_score, disc):
        #read the bioassays first
        assay1 = read_assay_feature(bioAssay1)
        assay2 = read_assay_feature(bioAssay2)

        #initialize parameters
        align_mode = "fea"

	#match/mismatch/gap scores are not used. Instead, use tanimoto similarity
	match_score	= None
	mismatch_score	= None
	#gap_score	= None

        #calculate similarity
        similarity = sequence_alignment_similarity(assay1, assay2, align_mode, match_score, mismatch_score, gap_score, disc)

        return similarity




#calculate similarity by cross ranking based sequence alignment by compound ID
#use bioAssay1's model to predict bioassay2 training file, 
#then align by bioassay2's training file and bioassay2's prediction with bioassay1's model.
#bioAssay1= path to bioAssay1
#bioAssay2= path to bioAssay2
#model_file= path to save bioassay 1's model
#pred_file= path to save bioassay 2' prediction with bioassay 1's model
#c= parameter '-c' for SVM light to train bioassay 1's model
#t= parameter '-t' for SVM light to train bioassay 1's model
#match_score= user-defined score for a matched step
#mismatch_score= user-defined score for a mismatched step
#gap_score= user-defined score for a gap step
#disc = (bool) if True, alignment score is discounted by position
def x_cid_sim(bioAssay1, bioAssay2, model_file, pred_file, c, t, match_score, mismatch_score, gap_score, disc):
        #read the bioassays first
        assay1 = read_assay_feature(bioAssay1)
        assay2 = read_assay_feature(bioAssay2)

	#call cross ranking based method to train bioAssay1's model
        #and predict bioAssay2 with bioAssay1's model
        cross_ranking_method(bioAssay1, bioAssay2, model_file, pred_file, c, t)

	#read ground truth scores and prediction scores
	true_score = read_assay_score(bioAssay2)
	pred_score = read_assay_score(pred_file)

	#sort scores from high to low
        rankT = sorted(range(len(true_score)), key=lambda k:true_score[k], reverse=True)
        rankP = sorted(range(len(pred_score)), key=lambda k:pred_score[k], reverse=True)

        #re-arrange features based on sorted scores
        sorted_assay1 = [assay1[i] for i in rankT]
        sorted_assay2 = [assay2[i] for i in rankP]


        #initialize parameters
        align_mode = "cid"

        #calculate similarity
        similarity = sequence_alignment_similarity(sorted_assay1, sorted_assay2, align_mode, match_score, mismatch_score, gap_score, disc)

        return similarity




#calculate similarity by cross ranking based sequence alignment by compound similarity
#use bioAssay1's model to predict bioassay2 training file, 
#then align by bioassay2's training file and bioassay2's prediction with bioassay1's model.
#bioAssay1= path to bioAssay1
#bioAssay2= path to bioAssay2
#model_file= path to save bioassay 1's model
#pred_file= path to save bioassay 2' prediction with bioassay 1's model
#c= parameter '-c' for SVM light to train bioassay 1's model
#t= parameter '-t' for SVM light to train bioassay 1's model
#disc = (bool) if True, alignment score is discounted by position
def x_csim_sim(bioAssay1, bioAssay2, model_file, pred_file, c, t, gap_score, disc):
        #read the bioassays first
        assay1 = read_assay_feature(bioAssay1)
        assay2 = read_assay_feature(bioAssay2)

        #call cross ranking based method to train bioAssay1's model
        #and predict bioAssay2 with bioAssay1's model
        cross_ranking_method(bioAssay1, bioAssay2, model_file, pred_file, c, t)

        #read ground truth scores and prediction scores
        true_score = read_assay_score(bioAssay2)
        pred_score = read_assay_score(pred_file)

        #sort scores from high to low
        rankT = sorted(range(len(true_score)), key=lambda k:true_score[k], reverse=True)
        rankP = sorted(range(len(pred_score)), key=lambda k:pred_score[k], reverse=True)

        #re-arrange features based on sorted scores
        sorted_assay1 = [assay1[i] for i in rankT]
        sorted_assay2 = [assay2[i] for i in rankP]


        #initialize parameters
        align_mode = "fea"
	
	#match/mismatch/gap scores are not used. Instead, use tanimoto similarity
        match_score     = None
        mismatch_score  = None
        #gap_score       = None


        #calculate similarity
        similarity = sequence_alignment_similarity(sorted_assay1, sorted_assay2, align_mode, match_score, mismatch_score, gap_score, disc)

        return similarity





#calcualte the similarity between two bio assays
def BA_sim(bioAssay1, bioAssay2, mode, model_file, pred_file, c, t, match_score, mismatch_score, gap_score):
	#initialize similarity
	similarity = 0.0
	
	#check which method to use
	#mode "pairwise": 
	#use average of compound pair-wise similarity
	if mode == 'pairwise':
		similarity = pairwise_sim(bioAssay1, bioAssay2)

	#mode "xci": 
	#use cross ranking based CI method 
	elif mode == 'xci':
		similarity = xci_sim(bioAssay1, bioAssay2, model_file, pred_file, c, t)

	#mode "p_cid_disc": 
	#profiling based sequence alignment with compound identities, with positional discount
	elif mode == 'p_cid_disc':
		similarity = p_cid_sim(bioAssay1, bioAssay2, match_score, mismatch_score, gap_score, "disc")

	#mode "p_cid_undisc": 
        #profiling based sequence alignment with compound identities
        elif mode == 'p_cid_undisc':
		similarity = p_cid_sim(bioAssay1, bioAssay2, match_score, mismatch_score, gap_score, "undisc")

	#mode "p_csim_disc": 
        #profiling based sequence alignment with compound similarities, with positional discount
        elif mode == 'p_csim_disc':
		similarity = p_csim_sim(bioAssay1, bioAssay2, gap_score, "disc")

        #mode "p_csim_undisc": 
        #profiling based sequence alignment with compound similarities
        elif mode == 'p_csim_undisc':
		similarity = p_csim_sim(bioAssay1, bioAssay2, gap_score, "undisc")

	#mode "x_cid_disc": 
        #cross ranking based sequence alignment with compound identities, with positional discount
        elif mode == 'x_cid_disc':
		similarity = x_cid_sim(bioAssay1, bioAssay2, model_file, pred_file, c, t, match_score, mismatch_score, gap_score, "disc")

        #mode "x_cid_undisc": 
        #cross ranking based sequence alignment with compound identities
        elif mode == 'x_cid_undisc':
		similarity = x_cid_sim(bioAssay1, bioAssay2, model_file, pred_file, c, t, match_score, mismatch_score, gap_score, "undisc")

        #mode "x_csim_disc": 
        #cross ranking based sequence alignment with compound similarities, with positional discount
        elif mode == 'x_csim_disc':
		similarity = x_csim_sim(bioAssay1, bioAssay2, model_file, pred_file, c, t, gap_score, "disc")

        #mode "x_csim_undisc": 
        #cross ranking based sequence alignment with compound similarities
        elif mode == 'x_csim_undisc':
		similarity = x_csim_sim(bioAssay1, bioAssay2, model_file, pred_file, c, t, gap_score, "undisc")

	#other method input are invalid
	else:
		sys.exit("Please input a supported method to calculate bioassay similarities.")

		
	return similarity



