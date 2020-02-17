#format from converted to training and test files
import csv
import os
import math
import sys


def CI(truth, pred):
	truthArr = []
	predArr  = []

	
	#construct an array of truth scores
	with open(truth) as t:
		#for each row, the ground truth score is the first column
		for row in t:
			rowt = row.split(' ')
			truthArr.append(float(rowt[0]))

	#construct an array of pred scores respective to the truth records
	with open(pred) as p:
		#for each row, the prediction score is the whole line
		for row in p:
			predArr.append(float(row.rstrip('\n')))
	
	#initialze some metrics
	trueCount = 0.0
	size  = len(truthArr)
	CI_val = 0.0
	
	#if size is less than 1, it is not meaningful to calculate a CI
	if size > 1:
		#count true ranking pairs
		for i in range (0, size - 1):
			for j in range(i + 1, size):
				#if rankings of two items are consistant in ground truth and prediction
				#then it is a correctly ranking pair
				if((truthArr[i] > truthArr[j]) and (predArr[i] > predArr[j])):
                                        trueCount += 1
                                elif((truthArr[i] < truthArr[j]) and (predArr[i] < predArr[j])):
                                        trueCount += 1
                                elif((truthArr[i] == truthArr[j]) and (predArr[i] == predArr[j])):
                                        trueCount += 1

		#calculate CI
		CI_val = float(trueCount)/(0.5*size*(size-1))
	
	else:
		sys.exit('Size too small. Unable to calculate CI.')	
	
	return CI_val
	
		



def write_CI(truth_file, pred_file, CI_file):
	
	#set parameters and call rankerror function
	CI_val = CI(truefile , predfile)
	#open file writer
	f_new = open(CI_file, 'w')
	#write value
	f_new.write('%6.4f\n' % (CI_val))
	#close file
	f_new.close()


