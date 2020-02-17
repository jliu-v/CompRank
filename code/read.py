import sys



#read an assay feature from assay file
#assay_file= path to assay file
def read_assay_feature(assay_file):
        #initialize empty array to store compounds
        features = []

        #open file reader
        with open(assay_file, 'r') as reader:
                #rows appended are one line in bioassay
                #online contains score, qid, and feature
                #line to be parsed in next step
                for row in reader:
                        parts = row.strip("\n").split()

                        #from the third column are the features
                        #only features are needed to calculate similarities
                        features.append(parts[2:])
        #close file reader
        reader.close()

        #return the stored assay
        return features


#read an assay scores from assay file
#assay_file= path to assay file
def read_assay(assay_file):
        #initialize empty array to store compounds
        scores  = []

        #open file reader
        with open(assay_file, 'r') as reader:
                #rows appended are one line in bioassay
                #online contains score, qid, and feature
                #line to be parsed in next step
                for row in reader:
                        parts = row.strip("\n").split()

                        #the scores are always in the first column
                        scores.append(float(parts[0]))

        #close file reader
        reader.close()

        #return the stored assay
        return scores



#read an assay from assay file
#assay_file= path to assay file
def read_assay(assay_file):
	#initialize empty array to store compounds
	features = []
	scores  = []

	#open file reader
	with open(assay_file, 'r') as reader:
                #rows appended are one line in bioassay
                #online contains score, qid, and feature
                #line to be parsed in next step
                for row in reader:
			parts = row.strip("\n").split()

			#the scores are always in the first column
			scores.append(float(parts[0]))

                        #from the third column are the features
                        #only features are needed to calculate similarities
                        features.append(parts[2:])
        #close file reader
	reader.close()
	
	#return the stored assay
	return {"features":features, "scores":scores}




