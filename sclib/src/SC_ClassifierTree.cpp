/**************************************************************************/
/*    Responsibility:																											*/
/*      - A binary tree of classifiers as an alternative way to do multi- */
/*        class classification as suggested in "Content-based audio       */
/*        classification and segmentation by using support vector         */
/*        machines", Lu, Zhang, Li, 2003                                  *
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.05.2006																								*/
/**************************************************************************/

#include <limits.h>
#include "SC_ClassifierTree.h"
#include "SC_Aux.h"
#include "SC_GroundTruth.h"
#include "SC_ClassifierHandler.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_ClassifierTree::SC_ClassifierTree(SC_TweakableParameters* pTweak, bool doScaling, bool verbose, bool traverseCompleteTree, int nodeCount, int **treeConstructionTable, int **posLabels, int *posLabelsCount, int **negLabels, int *negLabelsCount, int *posResult, int *negResult, unsigned int classifierType) : SC_Classifier(pTweak, doScaling) {
	this->classifierType = sclib::ctTree;
	this->classifierType = classifierType;
  this->verbose = verbose;
	this->traverseCompleteTree = traverseCompleteTree;

  if (nodeCount == 0) {
    this->pRoot = NULL;
		this->pFlatTree = NULL;
  } else {
    this->pRoot = constructTree(nodeCount, treeConstructionTable, posLabels, posLabelsCount, negLabels, negLabelsCount, posResult, negResult, classifierType);
  }
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_ClassifierTree::~SC_ClassifierTree() {
  MFree_0D(this->pRoot);
	MFree_1D(this->pFlatTree);
}

//====================================================================================================================
//	construct a binary tree (not necessarily balanced) out of the supplied data:
//
//  nodeCount: nr of nodes (including root and leaves) in the tree; all following arrays/matrixes should have as many 
//             rows as there are nodes; the indices into these structures are used as the nodeNr for the respective 
//             nodes
//  treeConstructionTable: nodeNr as index, array according to |nodeNr of pos-sibling|nodeNr of neg-sibling|; -1 as 
//                         nodeNr of a sibling indicates no sibling at this position. An example:
//                                    0             [0]->| 1| 2|
//                                  /   \      =>   [1]->|-1|-1|
//                                 1     2          [2]->| 3| 4|
//                                      / \         [3]->|-1|-1|
//                                     3   4        [4]->|-1|-1|
//                         siblings must have greater nodeNrs than their parents (greater as the index of the row in 
//                         which they appear); the root should have nodeNr/index 0
//  posLabels: list of labels used for the positive class for the node who's nodeNr corresponds with the index
//  posLabelsCount: length of that list
//  negLabels: s.a.
//  negLabelsCount: s.a.
//  posResults: the label that should be returned from that node if it is a leave and classifies a sample positively
//  negResult: s.a.
//  classifiyerType: an SCLIB_CT_* constant to choose between SVM, ML, ...
//====================================================================================================================
SC_ClassifierTree::SC_ClassifierTreeNode* SC_ClassifierTree::constructTree(int nodeCount, int **treeConstructionTable, int **posLabels, int *posLabelsCount, int **negLabels, int *negLabelsCount, int *posResult, int *negResult, unsigned int classifierType) {
  long int node;
  SC_ClassifierTree::SC_ClassifierTreeNode *pRoot = NULL;
  SC_Classifier *pClassifier = NULL;

  if (nodeCount > 0) {
		this->nodeCount = nodeCount;

    //create space
    MArray_1D(this->pFlatTree, nodeCount, SC_ClassifierTree::SC_ClassifierTreeNode*, "SC_ClassifierTree.constructTree: pFlatTree");

    for (node = nodeCount-1; node >= 0; node--) { //construct the leave first, went on to the root
      //build the node itself
      this->pFlatTree[node] = new SC_ClassifierTree::SC_ClassifierTreeNode(this->pTweak, node, this->classifierType, this->verbose, posLabels[node], posLabelsCount[node], negLabels[node], negLabelsCount[node], posResult[node], negResult[node]);

      //test the assumptions
      if ((treeConstructionTable[node][0] <= node && treeConstructionTable[node][0] != sclib::noNode) || (treeConstructionTable[node][1] <= node && treeConstructionTable[node][1] != sclib::noNode)) {
        REPORT_ERROR(SVLIB_BadArg, "Error in the tree-construction table: Siblings must have higher nodeNrs than their parents!");
      }

      //connect it's siblings (they have been previously created due to the "siblings have higher nodeNrs than their parents"-assumption
      if (treeConstructionTable[node][0] >= 0) { //left = positive
        this->pFlatTree[node]->setLeftNode(this->pFlatTree[treeConstructionTable[node][0]]);
      } else {
        this->pFlatTree[node]->setLeftNode(NULL);
      }
      if (treeConstructionTable[node][1] >= 0) { //right = negative
        this->pFlatTree[node]->setRightNode(this->pFlatTree[treeConstructionTable[node][1]]);
      } else {
        this->pFlatTree[node]->setRightNode(NULL);
      }
    }
    pRoot = this->pFlatTree[0];
  }

  return pRoot;
}

//====================================================================================================================
//	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
//====================================================================================================================
int SC_ClassifierTree::trainTwoClass(SV_Data *pPositive, SV_Data *pNegative) {
  int *classes, res = SVLIB_Fail;
  SV_Data *pCompleteData = NULL, *pHook = pPositive->Next;

  if (this->verbose == true) {
    printf("%s\n", "Complaint: Two-class training in case of the classification-tree makes no sense at all... doing it anyway.");
  }

  pPositive->Next = pNegative;
  pCompleteData = pPositive->MergeData(2);
  MArray_1D(classes, pCompleteData->Row, int, "SC_ClassifierTree.TrainTwoClasses: classes");
  for (long int x = 0; x < pCompleteData->Row; x++) {
    classes[x] = (x < pPositive->Row) ? sclib::labelPositive : sclib::labelNegative;
  }
  res = trainMultiClass(pCompleteData, classes);
  MFree_0D(pCompleteData);
  MFree_1D(classes);
  pPositive->Next = pHook;

  return res;
}

//====================================================================================================================
//	train a classifier for distinguishing between several classes
//  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
//  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
//  respective row of pData.
//====================================================================================================================
int SC_ClassifierTree::trainMultiClass(SV_Data *pData, int *classes) {
  int res = SVLIB_Fail;
  SV_Data *pScaledData = NULL;
	SC_ClassifierHandler handler(this->pTweak);

  //mark current tree as untrained and free old scaling parameters; 
	//it is assumed that a (untrained) tree already exists as constructed by constructTree()
	//because tree construction is not handled here, the (old) classifier can't be freed here prior to (new) training
  MFree_0D(this->pScale);
  //MFree_0D(this->pRoot);
	//MFree_1D(this->pFlatTree);
  this->isTrained = false;
	this->classCount = 0;

	//test the implicit assumption
	if (handler.canHandleMulticlassData(this->classifierType) != true) {
		return SVLIB_Fail;
	}

  if (this->pRoot != NULL) {
    if (this->doScaling == true) { //find scaling parameters
      this->pScale = findScalingParameters(pData);
    }

    pScaledData = scaleFeatures(pData, this->pScale, -1, true); //save some memory by just linking in case of no scaling
    res = this->pRoot->train(pData, classes);
    if (pScaledData != pData) {
      MFree_0D(pScaledData);
    }
  
    if (res == SVLIB_Fail) {
      this->isTrained = false;
      MFree_0D(this->pScale);
    } else {
      this->isTrained = true;
			this->classCount = getDistinctClassCount(classes, pData->Row);
    }
  }

  return res;
}

//====================================================================================================================
//	classifiy previously unseen test-data; returned is an array of classlabels, each entry corresponding to the 
//  respective row in pData; if available, the probabilities for each class-decision are given in the pProbabilities 
//  parameter, but here, the array-indexing has a different meaning than in the other classifier classes:
//  in col. nodeNr*2 is the label/labels assigned by this node, in nodeNr*2+1 is the correspondig probability
//
//  ATTENTION: only a previously trained (or loaded) classifier can be used for classification!
//====================================================================================================================
int* SC_ClassifierTree::classify(SV_Data *pData, SV_Data* &pProbabilities) {
  int *classes = NULL;
  SV_Data *pScaledData = NULL;
  
  if (this->doScaling == true && this->pScale == NULL) {
    REPORT_ERROR(SVLIB_BadData, "Can't do scaling if no scaling parameters where found");
  }

  if (this->isTrained == true) {
    if (this->pRoot != NULL) {
      
			//here, the array-indexing has a different meaning than in the other classifier classes:
			//in col. nodeNr*2 is the label assigned by this node, in nodeNr*2+1 is the correspondig probability
			MFree_0D(pProbabilities);
			pProbabilities = new SV_Data(pData->Row, this->nodeCount*2); 
      MArray_1D(classes, pData->Row, int, "SC_ClassifierTree.classify: classes");
      
      for (long int y = 0; y < pData->Row; y++) { //classify each feature vector separately
        pScaledData = scaleFeatures(pData, this->pScale, y);
        this->pRoot->setBothChildrenMode(this->traverseCompleteTree);
        classes[y] = this->pRoot->classify(pScaledData, pProbabilities->Mat[y]);
        MFree_0D(pScaledData);
      }

    } else {
      this->isTrained = false;
    }
  }

  return classes;
}

//====================================================================================================================
//	save a trained classifier to a file
//====================================================================================================================
int SC_ClassifierTree::saveClassifier(const char *fileName) {
  int res = SVLIB_Fail;
  char *scaleFileName = NULL;
  SV_DataIO io;

  if (strlen(fileName) > 0) {
    if (this->isTrained == true) {
      if (this->pRoot != NULL) {
        res = this->pRoot->save(fileName);
        if (res == SVLIB_Ok && this->doScaling == true) { //save also the scaling parameters of this training set to use it with the test-data to come
          scaleFileName = sclib::exchangeFileExtension(fileName, ".scale");
          io.OpenFile(scaleFileName, WRITE_REC);
		      res = io.PutDataRec(*this->pScale);
          res = (res == 0) ? SVLIB_Fail : SVLIB_Ok; //translate between different traditions to report errrors or success...
          io.CloseFile();
          MFree_1D(scaleFileName);
        }
      } else {
        this->isTrained = false;
      }
    }
  }

  return res;
}

//====================================================================================================================
//	load a trained classifier from a file
//====================================================================================================================
int SC_ClassifierTree::loadClassifier(const char *fileName) {
  int res;
  char *scaleFileName = NULL;
  SV_DataIO io;
	int *classes = NULL, counter = 0;

  //load classifier
	this->nodeCount = 0;
  MFree_0D(this->pRoot);
	MFree_1D(this->pFlatTree);
  this->pRoot = new SC_ClassifierTree::SC_ClassifierTreeNode(this->pTweak, 0, this->classifierType); //create the raw root node
  res = this->pRoot->load(fileName); //this loads the root and all it's siblings
  if (this->pRoot != NULL && res != SVLIB_Fail) {
		this->nodeCount = getTreeSize(this->pRoot);
		MArray_1D(this->pFlatTree, this->nodeCount, SC_ClassifierTree::SC_ClassifierTreeNode*, "SC_ClassifierTree.loadClassifier: pFlatTree");
		fillFlatTree(this->pRoot, this->pFlatTree);
    res = SVLIB_Ok;
    this->isTrained = true;
  } else {
    res = SVLIB_Fail;
    this->isTrained = false;
  }

  //also try to load (if exists) the scaling parameters:
  MFree_0D(this->pScale);
  if (res == SVLIB_Ok && this->doScaling == true) {
    scaleFileName = sclib::exchangeFileExtension(fileName, ".scale");
    if (sclib::fileExists(scaleFileName) == true) {
      io.OpenFile(scaleFileName, READ_REC);
      this->pScale = io.GetAllRec();
      io.CloseFile();
    }
    MFree_1D(scaleFileName);
  }

  //correct a bug that may have occured during (very time-consuming, so retraining is a bad alternatice) classifier-training in earlier versions: 
  //the assignment of positive and negative nodes to left and right siblings was different inside the algorithms and on the caller side, 
  //so each node is now correctly trained, but the tree is wrongly connected
  if (this->isTrained == true) {
    if (resultContradiction(this->pRoot) == true) { //try to guess is the tree is wrongly connected
      printf("The loaded classifier seems to be wrongly connected due to an old bug -> corrected by switching pos/neg nodes.\n");
      switchLeftRightNodes(this->pRoot);
    }
  }

	//fill the classCount-member by deducing nr. of classes from tree-leaves
	MArray_1D(classes, this->nodeCount*2, int, "SC_ClassifierTree:loadClassifier: classes");
	for (int n = 0; n < this->nodeCount; n++) {
		if (this->pFlatTree[n]->getLeftNode() == NULL && this->pFlatTree[n]->getRightNode() == NULL) { //a leaf...
			classes[counter++] = this->pFlatTree[n]->getPositiveResult();
			classes[counter++] = this->pFlatTree[n]->getNegativeResult();
		}
	}
	this->classCount = getDistinctClassCount(classes, counter);
	MFree_1D(classes);

  return res;
}

//====================================================================================================================
//	to convert between given labels and indices into the probability-parameter of the classifiy()-method
//  returns -1 if something goes wrong
//====================================================================================================================
long int SC_ClassifierTree::label2idx(long int label) {
	long int res = -1;

	return res; //returns error, cause the label is better deduced from the pProbabilities-object itself in this class
}

//====================================================================================================================
//	to convert between given labels and indices into the probability-parameter of the classifiy()-method
//  returns sclib::noType if something goes wrong
//====================================================================================================================
long int SC_ClassifierTree::idx2label(long int idx) {
	long int res = sclib::noType;
	
	return res;  //returns error, cause the idx is better deduced from the pProbabilities-object itself in this class
}

//====================================================================================================================
//	traverses the tree recursively and returns the overall nodeCount
//====================================================================================================================
int SC_ClassifierTree::getTreeSize(SC_ClassifierTree::SC_ClassifierTreeNode *pStart) {
	int leftCount, rightCount;

	leftCount = (pStart->getLeftNode() != NULL) ? getTreeSize(pStart->getLeftNode()) : 0;
	rightCount = (pStart->getRightNode() != NULL) ? getTreeSize(pStart->getRightNode()): 0;

	return 1 + leftCount + rightCount;
}

//====================================================================================================================
//	traverses the tree recursively and successively fills the pointers in pFlatTree, which has to be allocated 
//  beforehand
//====================================================================================================================
void SC_ClassifierTree::fillFlatTree(SC_ClassifierTree::SC_ClassifierTreeNode *pStart, SC_ClassifierTree::SC_ClassifierTreeNode **pFlatTree) {
  pFlatTree[pStart->getNodeNr()] = pStart;

	if (pStart->getLeftNode() != NULL) {
		fillFlatTree(pStart->getLeftNode(), pFlatTree);
	}

	if (pStart->getRightNode() != NULL) {
		fillFlatTree(pStart->getRightNode(), pFlatTree);
	}
	
	return;
}

//====================================================================================================================
//	this method corrects in error that may have been occured during training: if in the tree-construction table 
//  positive- and negative nodes where switched, this can be corrected with this method for an already trained, but
//  wrongly connected classifier-tree
//====================================================================================================================
void SC_ClassifierTree::switchLeftRightNodes(SC_ClassifierTree::SC_ClassifierTreeNode *pNode) {
  SC_ClassifierTree::SC_ClassifierTreeNode *pHook = this->pRoot, *pTemp = NULL;

  if (pHook != NULL) {
    //exchange the two connections with each other
    pTemp = pHook->getLeftNode();
    pHook->setLeftNode(pHook->getRightNode(), false);
    pHook->setRightNode(pTemp, false);

    //recursively traverse the tree
    switchLeftRightNodes(pHook->getLeftNode());
    switchLeftRightNodes(pHook->getRightNode());
  }
  
  return;
}

//====================================================================================================================
//	returns true if inconsistencies seem to exist in the connection of the nodes according to the results of parents
//  and theire siblings
//  TODO: assumes result-types to be bitflags in 32bit variables
//====================================================================================================================
bool SC_ClassifierTree::resultContradiction(SC_ClassifierTree::SC_ClassifierTreeNode *pNode) {
  bool contradicts = false;
  long int x, y, rootLabel, rootRes, sibLabel, sibRes;
  
  //assume switched siblings if the positive result of the root is diametral to the positive result of the positive (left) sibling
  if (pNode->getLeftNode() != NULL) {
    rootRes = pNode->getPositiveResult();
    sibRes = pNode->getLeftNode()->getPositiveResult();
    rootLabel = 1; 
    for (x = 0; x < 32; x++) { //decompose result-types (may be combined ones, we need single bitflags for the isOppositeType()-test)
      sibLabel = 1;
      if (sclib::bitTest(rootRes, rootLabel) != 0) { //the type rootLabel is present in the positive result of the root
        for (y = 0; y < 32; y++) {
          if (sclib::bitTest(sibRes, sibLabel) != 0) { //the type sibLabel is present in the positive result of the left sibling
            if (SC_GroundTruth::isOppositeType(sibLabel, rootLabel) == true) {
              contradicts = true;
              break;
            }
          }
          sibLabel *= 2;
        }
      }
      if (contradicts == true) {
        break;
      }
      rootLabel *= 2;
    }
    if (contradicts == false) { //traverse the tree deeper if no contradiction found till here (depth first search left-2-right)
      contradicts = resultContradiction(pNode->getLeftNode());
    }
  }

  if (contradicts == false && pNode->getRightNode() != NULL) { //everything seems ok with the positive (left) branch, so check the right/negative one
    rootRes = pNode->getNegativeResult();
    sibRes = pNode->getRightNode()->getNegativeResult();
    rootLabel = 1; 
    for (x = 0; x < 32; x++) { //decompose result-types (may be combined ones, we need single bitflags for the isOppositeType()-test)
      sibLabel = 1;
      if (sclib::bitTest(rootRes, rootLabel) != 0) { //the type rootLabel is present in the negative result of the root
        for (y = 0; y < 32; y++) {
          if (sclib::bitTest(sibRes, sibLabel) != 0) { //the type sibLabel is present in the negative result of the right sibling
            if (SC_GroundTruth::isOppositeType(sibLabel, rootLabel) == true) {
              contradicts = true;
              break;
            }
          }
          sibLabel *= 2;
        }
      }
      if (contradicts == true) {
        break;
      }
      rootLabel *= 2;
    }
    if (contradicts == false) { //traverse the tree deeper if no contradiction found till here (depth first search left-2-right)
      contradicts = resultContradiction(pNode->getRightNode());
    }
 }

  return contradicts;
}
