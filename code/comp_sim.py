
#input compound feature vector, calculatet the similarity between two compounds
def tanimoto(arr1, arr2):
	size1      = float(len(arr1))
	size2      = float(len(arr2))
	sizeCommon = float(len(set(arr1).intersection(arr2)))
	

	if size1+size2-sizeCommon == 0:
		return 0
	return (sizeCommon/(size1+size2-sizeCommon))





