import numpy, ase

def writexyzatoms(atoms,filename, w='w',comment='Produced in TSASE'):
	xyz = open(filename,w)
	positions = atoms.get_positions()
	symbols = atoms.get_chemical_symbols()
	energy = atoms.get_potential_energy()
	xyz.write(str(len(positions)) + " \n")
	xyz.write(comment + " \n")
	for i in range(len(positions)):
		xyzposition = symbols[i] + " " + str(positions[i,0]) + " " + str(position[i,1]) + " " + str(position[i,2]) + " \n" 
		xyz.write(xyzposition)
	xyz.close()
	
def readxyzatoms(filename):
	xyz = open(filename, 'r') 
	numatoms = int(xyz.readline())
	dataY = [float(xyz.readline().split()[0])]
	datapoints = []
	for i in range(numatoms):
		vals = xyz.readline().split()
		vals.pop(0)
		vals[0] = float(vals[0])
		vals[1] = float(vals[1])
		vals[2] = float(vals[2])
		datapoints.append(vals)
	dataX = [numpy.asarray(datapoints).ravel()]
	line = xyz.readline().split() # check for EOF
	while (line != []): # this will read until EOF 	
		dataY.append(float(xyz.readline().split()[0]))
		datapoints = []
		for i in range(numatoms):
			vals = xyz.readline().split() 
			vals.pop(0)
			vals[0] = float(vals[0])
			vals[1] = float(vals[1])
			vals[2] = float(vals[2])
			datapoints.append(vals)
		dataX.append(numpy.asarray(datapoints).ravel())
		line = xyz.readline().split() # check for EOF
	dataXa = numpy.array(dataX)
	dataYa = numpy.array(dataY)
	return dataXa, dataYa
	
def learn(parameters,datapoints,datalabels,datapointstest,datalabelstest,kval=2): # Returns selSVM
	selSVM = selectmodel(parameters,datapoints,datalabels,kval)
	outValue = selSVM.decision_function(datapointstest) # don't forget that outValue is for the test set!
	out = selSVM.predict(datapointstest)
	try:
		from sklearn import metrics as me
	except:
		from scikits.learn import metrics as me
	print me.classification_report(datalabelstest,out)
	return selSVM

def selectmodel(parameters,xtrain,ytrain,kval):
	try:
		from sklearn import svm
		from sklearn.grid_search import GridSearchCV 
		from sklearn import metrics as me
		from sklearn.cross_validation import StratifiedKFold
	except:
		from scikits.learn import svm
		from scikits.learn.grid_search import GridSearchCV 
		from scikits.learn import metrics as me
		from scikits.learn.cross_val import StratifiedKFold
	from pprint import pprint
 
	mySVM = svm.SVC() 
	clf = GridSearchCV(mySVM, parameters, score_func=me.zero_one_score, 
	verbose=0, cv = StratifiedKFold(ytrain, k=kval)) 
	out = clf.fit(xtrain, ytrain).predict(xtrain)
	bestSVM = clf.best_estimator

	print "Report for the best estimator: " 
	print clf.best_estimator
	pprint(clf.grid_scores_)
	print
	return bestSVM

def randomize_datapoints(datapoints,datalabels):
	length = len(datapoints)
	aslist = datapoints.tolist()
	aslistlabel = datalabels.tolist()
	for i in range(length * 4):
		rando = numpy.random.randint(length)
		point = aslist.pop(rando)
		label = aslistlabel.pop(rando)
		aslist.append(point)
		aslistlabel.append(label)
	newdata = numpy.asarray(aslist)
	newlabels = numpy.asarray(aslistlabel)
	return newdata, newlabels
