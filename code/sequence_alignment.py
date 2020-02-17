from __future__ import print_function
from math import *
from comp_sim import *

def match(arr1, arr2, mode, match_score, mismatch_score):
	
	if   mode == 'cid':
		if(arr1 == arr2):

			return match_score
		else:
			return mismatch_score
	
	
	elif mode == 'fea':
		return tanimoto(arr1, arr2)
		
		
def sequence_alignment_similarity(seq1, seq2, mode, match_score, mismatch_score, gap_score, disc):
	
	#matrix to save datapath and result
	m   = len(seq1) + 1
	n   = len(seq2) + 1
	SIM = [[ 0 for col in range(n)] for row in range(m)]
	MAP = [['' for col in range(n)] for row in range(m)]
	
	#initialize worst cast
	for i in range(0, m):
		SIM[i][0] = i * gap_score
		MAP[i][0] = 'U'
	
	for j in range(1, n):
		SIM[0][j] = j * gap_score
		MAP[0][j] = 'L'
	
	#initialize positional discounts
	#if not undiscounted, discount will be 1	
	discount_i = 1.0
	discount_j = 1.0

	#updating worst case with dynamic programming
	for i in range(1, m):
		for j in range(1, n):
			#choose a max score among three options
			#seq1 gap, seq2 gap, or no gap
			
			if disc == 'disc':
				discount_i = 0.5 + (1 / (1 + exp( 0.1*(i-1) ) ) ) 
				discount_j = 0.5 + (1 / (1 + exp( 0.1*(j-1) ) ) )
			
			#U: update from "upside" to "downside"
			U  = SIM[i-1][j] + discount_i * gap_score
			#L: update from "left" to "right"
			L  = SIM[i][j-1] + discount_j * gap_score
			#UL: update from "upper left" side
			UL = SIM[i-1][j-1] + sqrt( (discount_i*discount_j) ) * match(seq1[i-1], seq2[j-1], mode, match_score, mismatch_score)

	 

			#keep track of path
			#which direction is the previous step
			max_num  = U
			max_map  = 'U'
			
			if(max_num < L):
				max_num = L
				max_map = 'L'
				
			if(max_num < UL):
				max_num = UL
				max_map = 'UL'
			
			#update the element in score matrix and step matrix
			SIM[i][j] = max_num
			MAP[i][j] = max_map
	
	
	#track back the optimal path
	#commented out because this part is not used in the project protocol
	#but may be useful in degub or future development
	'''
	#trace back for to generate back-path
	path1 = []
	path2 = []
	
	i = m - 1
	j = n - 1
	

	#keep going until hit the origin
	while(i > 0 or j > 0):
		#get current location
		current_position  = MAP[i][j]
		current_step_1    = []
		current_step_2    = []
		
		choice = 0
		
		#choose one from the three ways
		#gap on seq2
		if i >= 1 and current_position == 'U':
			current_step_1 = [i]
			current_step_2 = ['-']
			choice = 1
		
		#gap on seq1
		if j >= 1 and current_position == 'L':
			current_step_1 = ['-']
			current_step_2 = [j]
			choice = 2
		
		#no gap
		if i >= 1 and j >= 1 and current_position == 'UL':
			current_step_1 = [i]
			current_step_2 = [j]
			choice = 3
		
		#update path
		path1 = current_step_1 + path1
		path2 = current_step_2 + path2
		
		#update iterating indices
		if   choice == 1:
			i -= 1
		elif choice == 2:
			j -= 1
		elif choice == 3:
			i -= 1
			j -= 1
		else:
			#avoid infinite loop
			exit()
	#end while		
	
	#return the paths and the final similarity between the two sequence
	return [path1, path2, SIM[m-1][n-1]]
	'''
	
	#return the optimal alignment score
	return SIM[m-1][n-1]















