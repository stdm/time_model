/**************************************************************************/
/*    Responsibility:																											*/
/*      - A binary tree of classifiers as an alternative way to do multi- */
/*        class classification as suggested in "Content-based audio       */
/*        classification and segmentation by using support vector         */
/*        machines", Lu, Zhang, Li, 2003                                  */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.05.2006																								*/
/**************************************************************************/

#ifndef __SC_ClassifierTree_H__
#define __SC_ClassifierTree_H__

#include <iostream>
#include "SC_ClassifierHandler.h"
#include "SC_Classifier.h"
#include "SC_Api.h"
#include "SC_Aux.h"
#include "SC_TweakableParameters.h"
#include <SV_Data.h>
#include <SV_DataIO.h>

class SCLIB_API SC_ClassifierTree : public SC_Classifier {
  
  private :

    //====================================================================================================================
    //	auxiliary class to build the tree and work on it recursively
    //====================================================================================================================
    class SC_ClassifierTreeNode {
      private:
      protected:

        SC_Classifier *pClassifier; //one node in the tree itself

				//is it strange at first sight to assign the left node the positive label?!? imagine (computer-science) tree with the root as the top-node and the leaves below; it is quite natural to assign '+' to the left branch, then...
        SC_ClassifierTreeNode *pLeftNode; //this node further processes positively classified samples; NULL if the classifier is a leaf in the tree
        SC_ClassifierTreeNode *pRightNode; //s.a., negative=positive

        SC_TweakableParameters *pTweak;
        int *positiveLabels; //which labels in the training set should be used for positive examples
        int positiveLabelsCount; //how many different labels as positive examples are there?
        int *negativeLabels; //s.a., positive=negative
        int negativeLabelsCount; //s.a., positive=negative
        int positiveResult; //if this node is a leaf, this label will be assigned to a positively classified sample
        int negativeResult; //s.a., positive=negative
        int nodeNr; //a number to identify the node
        bool verbose; //to control console output
        bool evaluateBothChildren; //if true, both branches will be evaluated (not only the best path) to populate the probabilities-array
        unsigned int classifierType; //which classifier-type to use?

        //====================================================================================================================
        //  decide which classifier-type to build, then build a raw (untrained) classifier of that type
        //====================================================================================================================
        SC_Classifier* createClassifier(void) {
					SC_ClassifierHandler *pHandler = new SC_ClassifierHandler(this->pTweak);
          SC_Classifier *pClassifier = pHandler->buildClassifier(this->classifierType, this->pTweak, false, this->verbose); //scaling is done in the meta-classifier (the tree), not in the nodes

					MFree_0D(pHandler);

          return pClassifier;
        }

        //====================================================================================================================
        //  copy together just those samples in pData that have a label also listed in wantedLabels
        //====================================================================================================================
        SV_Data* filterData(SV_Data *pData, int *labels, int *wantedLabels, int wantedLabelsCount) {
          long int x, y, count = 0;
          int z;
          SV_Data *pFilteredData = NULL;
          
          //determine how many samples will be in the filtered result
          for (x = 0; x < pData->Row; x++) {
            for (z = 0; z < wantedLabelsCount; z++) {
              if (labels[x] == wantedLabels[z]) {
                count++;
                break;
              }
            }
          }

          //alloc memory
          pFilteredData = new SV_Data(count, pData->Col);

          //copy the relevant samples
          count = 0;
          for (x = 0; x < pData->Row; x++) {
            for (z = 0; z < wantedLabelsCount; z++) {
              if (labels[x] == wantedLabels[z]) {
                for (y = 0; y < pData->Col; y++) {
                  pFilteredData->Mat[count][y] = pData->Mat[x][y];
                }
                count++;
                break;
              }
            }
          }

          return pFilteredData;
        }

      public:

        //====================================================================================================================
        //  the constructor copys the lists, but only links to the other
        //====================================================================================================================
        SC_ClassifierTreeNode(SC_TweakableParameters *pTweak, int nodeNr, unsigned int classifierType, bool verbose = true, int *posLabels = NULL, int posLabelsCount = 0, int *negLabels = NULL, int negLabelsCount = 0, int posResult = 0, int negResult = 0, SC_ClassifierTreeNode *pLeftNode = NULL, SC_ClassifierTreeNode *pRightNode = NULL, bool evaluateBothChildren = false) {
          int i;

          this->positiveResult = posResult;
          this->positiveLabelsCount = posLabelsCount;
          this->positiveLabels = NULL;
          if (this->positiveLabelsCount > 0) {
            MArray_1D(this->positiveLabels, this->positiveLabelsCount, int, "SC_ClassifierTree::SC_ClassifierTreeNode: positiveLabels");
            for (i = 0; i < this->positiveLabelsCount; i++) {
              this->positiveLabels[i] = posLabels[i];
            }
          }

          this->negativeResult = negResult;
          this->negativeLabelsCount = negLabelsCount;
          this->negativeLabels =  NULL;
          if (this->negativeLabelsCount > 0) {
            MArray_1D(this->negativeLabels, this->negativeLabelsCount, int, "SC_ClassifierTree::SC_ClassifierTreeNode: negativeLabels");
            for (i = 0; i < this->negativeLabelsCount; i++) {
              this->negativeLabels[i] = negLabels[i];
            }
          }

          this->pTweak = pTweak;
          this->classifierType = classifierType;
          this->verbose = verbose;
          this->nodeNr = nodeNr;
          this->pClassifier = createClassifier();
          this->pLeftNode = pLeftNode;
          this->pRightNode = pRightNode;
          this->evaluateBothChildren = evaluateBothChildren;
        }

        //====================================================================================================================
        //  the destructor recursively deletes all child-nodes and also the linked classifyer
        //====================================================================================================================
        virtual ~SC_ClassifierTreeNode() { 
          MFree_0D(this->pClassifier);
          MFree_1D(this->positiveLabels);
          MFree_1D(this->negativeLabels);

          MFree_0D(this->pLeftNode);
          MFree_0D(this->pRightNode);
        }

        //====================================================================================================================
        //  geter/seter
        //====================================================================================================================
        int getNodeNr(void) {return this->nodeNr;}
        bool getVerboseMode(void) {return this->verbose;}
        unsigned int getClassifierType(void) {return this->classifierType;}
        SC_ClassifierTreeNode* getRightNode(void) {return this->pRightNode;}
        SC_ClassifierTreeNode* getLeftNode(void) {return this->pLeftNode;}
        void setRightNode(SC_ClassifierTreeNode *pNode, bool freeOldNode = true) {
          if (freeOldNode == true) {
            MFree_0D(this->pRightNode)
          }; 
          this->pRightNode = pNode; 
          return;
        }
        void setLeftNode(SC_ClassifierTreeNode *pNode, bool freeOldNode = true) {
          if (freeOldNode == true) {
            MFree_0D(this->pLeftNode)
          }; 
          this->pLeftNode = pNode; 
          return;
        }
        int getPositiveResult(void) {return this->positiveResult;}
        int getNegativeResult(void){return this->negativeResult;}
        void setBothChildrenMode(bool evaluateBoth) {this->evaluateBothChildren = evaluateBoth; return;}

        //====================================================================================================================
        //  train the current node and all it's children
        //====================================================================================================================
        int train(SV_Data *pData, int *classes) {
					int res;
          double posRes = 0.0, negRes = 0.0; 
          SV_Data *pPositiveData = filterData(pData, classes, this->positiveLabels, this->positiveLabelsCount);
          SV_Data *pNegativeData = filterData(pData, classes, this->negativeLabels, this->negativeLabelsCount);

          res = this->pClassifier->trainTwoClass(pPositiveData, pNegativeData); //train this node
          MFree_0D(pPositiveData); //release copys of training data
          MFree_0D(pNegativeData);

          if (res == SVLIB_Fail) {
            REPORT_ERROR(SVLIB_Fail, "Error during node training");
          } else {
            if (this->pRightNode != NULL) { //recursively train right (positive) nodes
              res = this->pRightNode->train(pData, classes);
            }
            if (res != SVLIB_Fail && this->pLeftNode != NULL) { //recursively train left (negative) nodes
              res = this->pLeftNode->train(pData, classes);
            }
          }

          return res;
        }

        //====================================================================================================================
        //  classify the given sample (ATTENTION: there should be only one row in pSample!!!);
        //  call the child-nodes recursively for further classification if any, otherwise (this is a leave-node then) return 
        //  result directly
				//  in allLabelsAndProbs a float-array is given; each node in the tree has to save it's classification result (label)
				//  in the cell with index nodeNr*2 and the corresponding probability in nodeNr*2+1
				//  if evaluateBothChilds==true, not only the path following the predicted label is evaluated but both childnodes, in 
				//  order to get all probabilities listed in allLabelsAndProbs
        //====================================================================================================================
        int classify(SV_Data *pSample, float *allLabelsAndProbs) {
          int *result, res;
					SV_Data *pProb = NULL;
            
          if (pSample->Row != 1) {
            REPORT_ERROR(SVLIB_BadArg, "Each node can classify only one single sample at once");
          }

          result = this->pClassifier->classify(pSample, pProb);
					allLabelsAndProbs[this->nodeNr*2] = (float)(result[0]); //return only for each node if the result was positive or negative, not the label attached to this decision; this can be handled more flexible outside from the caller's side...
					allLabelsAndProbs[this->nodeNr*2 + 1] = pProb->Mat[0][this->pClassifier->label2idx(result[0])]; //probability of this label
					MFree_0D(pProb);

          if (result[0] == sclib::labelPositive) {
						MFree_1D(result);
            if (this->pLeftNode != NULL) {
              this->pLeftNode->setBothChildrenMode(this->evaluateBothChildren);
              res = this->pLeftNode->classify(pSample, allLabelsAndProbs);
            } else {
              res = this->positiveResult;
            }
						if (this->evaluateBothChildren == true && this->pRightNode != NULL) {
              this->pRightNode->setBothChildrenMode(this->evaluateBothChildren);
							this->pRightNode->classify(pSample, allLabelsAndProbs);
						}
          } else if (result[0] == sclib::labelNegative) {
						MFree_1D(result);
            if (this->pRightNode != NULL) {
              this->pRightNode->setBothChildrenMode(this->evaluateBothChildren);
              res = this->pRightNode->classify(pSample, allLabelsAndProbs);
            } else {
              res = this->negativeResult;
            }
						if (this->evaluateBothChildren == true && this->pLeftNode != NULL) {
              this->pLeftNode->setBothChildrenMode(this->evaluateBothChildren);
							this->pLeftNode->classify(pSample, allLabelsAndProbs);
						}
          } else {
            MFree_1D(result);
            REPORT_ERROR(SVLIB_BadData, "Classification-result is neither positive nor negative... check the result-types!");
          }

          return res;
        }
        
        //====================================================================================================================
        //  save a node (and all it's children) to a file; replace the file-extension with the string "node_" and the 
				//  respective nodeNr
        //====================================================================================================================
        int save(const char *fileName) {
          char *newFileName, newExtension[sclib::bufferSize];
          int bytes, res, leftNode = (this->pLeftNode != NULL) ? this->pLeftNode->getNodeNr() : sclib::noNode, rightNode = (this->pRightNode != NULL) ? this->pRightNode->getNodeNr() : sclib::noNode;
          fstream file;
					SV_DataIO io;
					SV_DataIO::SV_DatatypeSizes codeSizes;
					io.getCurrentDatatypeSizes(codeSizes);
          
          //save class information
          sprintf(newExtension, "%s%i\0", ".node_", this->nodeNr);
					newFileName = sclib::exchangeFileExtension(fileName, newExtension);
          file.open(newFileName, ios::out|ios::binary);
					bytes = io.writeMachineHeader(&file, codeSizes);
					bytes += io.writeScalar(&file, this->nodeNr);
					bytes += io.writeScalar(&file, this->verbose);
					bytes += io.writeScalar(&file, this->classifierType);
					bytes += io.writeScalar(&file, this->positiveLabelsCount);
					bytes += io.writeArray(&file, this->positiveLabels, this->positiveLabelsCount);
					bytes += io.writeScalar(&file, this->negativeLabelsCount);
					bytes += io.writeArray(&file, this->negativeLabels, this->negativeLabelsCount);
					bytes += io.writeScalar(&file, this->positiveResult);
					bytes += io.writeScalar(&file, this->negativeResult);
					bytes += io.writeScalar(&file, rightNode);
					bytes += io.writeScalar(&file, leftNode);
          res = (file.good() != TRUE) ? SVLIB_Fail : bytes;
	        file.close();

          if (res != SVLIB_Fail) {
            //save the classifier
            MFree_1D(newFileName);
            sprintf(newExtension, "%s%i\0", ".classifier_", this->nodeNr);
						newFileName = sclib::exchangeFileExtension(fileName, newExtension);
            res = this->pClassifier->saveClassifier(newFileName);
            MFree_0D(newFileName);

            if (res != SVLIB_Fail) {
              //save children
              if (this->pRightNode != NULL) {
                res = this->pRightNode->save(fileName);
              }
              if (res != SVLIB_Fail) {
                if (this->pLeftNode != NULL) {
                  res = this->pLeftNode->save(fileName);
                }
              }
            }
          }
          
          return (res != SVLIB_Fail) ? SVLIB_Ok : SVLIB_Fail;
        }

        //====================================================================================================================
        //  load a node (and all it's children) from a file; 
        //  it is assumed that it resides in a file with the forepart as given and the extension like the string "node_" and 
        //  the respective nodeNr the classifier is assumed to be in the file extended by "classifier_" & nodeNr
        //====================================================================================================================
        int load(const char *fileName) {
          char *newFileName, newExtension[sclib::bufferSize];
          int bytes, res, leftNode, rightNode;
          fstream file;
					SV_DataIO io;
					SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
					io.getCurrentDatatypeSizes(codeSizes);
          
          //destruct old information prior to loading
          MFree_1D(this->positiveLabels);
          MFree_1D(this->negativeLabels);
          MFree_0D(this->pLeftNode);
          MFree_0D(this->pRightNode);

          //load class information
          sprintf(newExtension, "%s%i\0", ".node_", this->nodeNr);
          newFileName = sclib::exchangeFileExtension(fileName, newExtension);
          file.open(newFileName, ios::in|ios::binary);
					bytes = io.readMachineHeader(&file, fileSizes, true);
					if (bytes > 0) {
						io.consumeBytes(&file, bytes);
					} else {
						bytes = 0;
					}
					bytes += io.readScalar(&file, this->nodeNr, codeSizes, fileSizes);
					bytes += io.readScalar(&file, this->verbose, codeSizes, fileSizes);
					bytes += io.readScalar(&file, this->classifierType, codeSizes, fileSizes);
					bytes += io.readScalar(&file, this->positiveLabelsCount, codeSizes, fileSizes);
          MArray_1D(this->positiveLabels, this->positiveLabelsCount, int, "SC_ClassifierTree::SC_ClassifierTreeNode.load: positiveLabels");
					bytes += io.readArray(&file, this->positiveLabels, this->positiveLabelsCount, codeSizes, fileSizes);
					bytes += io.readScalar(&file, this->negativeLabelsCount, codeSizes, fileSizes);
          MArray_1D(this->negativeLabels, this->negativeLabelsCount, int, "SC_ClassifierTree::SC_ClassifierTreeNode.load: negativeLabels");
					bytes += io.readArray(&file, this->negativeLabels, this->negativeLabelsCount, codeSizes, fileSizes);
					bytes += io.readScalar(&file, this->positiveResult, codeSizes, fileSizes);
					bytes += io.readScalar(&file, this->negativeResult, codeSizes, fileSizes);
					bytes += io.readScalar(&file, rightNode, codeSizes, fileSizes); //just the nodeNrs to know which file to load next
					bytes += io.readScalar(&file, leftNode, codeSizes, fileSizes);
          res = (file.good() != TRUE) ? SVLIB_Fail : bytes;
	        file.close();

          if (res != SVLIB_Fail) {
            //load the classifier
            MFree_1D(newFileName);
            sprintf(newExtension, "%s%i\0", ".classifier_", this->nodeNr);
            newFileName = sclib::exchangeFileExtension(fileName, newExtension);
            res = this->pClassifier->loadClassifier(newFileName);
            MFree_0D(newFileName);

            if (res != SVLIB_Fail) {
              //load children
              if (rightNode != sclib::noNode) {
                this->pRightNode = new SC_ClassifierTreeNode(this->pTweak, rightNode, this->classifierType);
                res = this->pRightNode->load(fileName);
              }
              if (res != SVLIB_Fail) {
                if (leftNode != sclib::noNode) {
                  this->pLeftNode = new SC_ClassifierTreeNode(this->pTweak, leftNode, this->classifierType);
                  res = this->pLeftNode->load(fileName);
                }
              }
            }
          }

          return (res != SVLIB_Fail) ? SVLIB_Ok : SVLIB_Fail;
        }
    }; //end of declaration/definition of SC_ClassifierTreeNode

  protected :
    
    SC_ClassifierTree::SC_ClassifierTreeNode *pRoot; //the root of the classification tree
		SC_ClassifierTree::SC_ClassifierTreeNode **pFlatTree; //to save pointers to each node directly, accessible through its nodeNr
    unsigned int classifierType;
    bool verbose; //to control console output...
		int nodeCount; //how many nodes are in the tree?
		bool traverseCompleteTree; //if true, not onyl the correct path along the tree is calculated, but all nodes are evaluated in order to get all probability results

    //====================================================================================================================
    //	auxiliary methods used during classifier loading
    //====================================================================================================================
		int getTreeSize(SC_ClassifierTree::SC_ClassifierTreeNode *pStart);
		void fillFlatTree(SC_ClassifierTree::SC_ClassifierTreeNode *pStart, SC_ClassifierTree::SC_ClassifierTreeNode **pFlatTree);

    //====================================================================================================================
    //	construct a binary tree (not necessarily balanced) out of the supplied data:
    //
    //  nodeCount: nr of nodes (including root and leaves) in the tree; all following arrays/matrixes should have as many 
    //             rows as there are nodes; the indices into these structures are used as the nodeNr for th respective 
    //             nodes
    //  treeConstructionTable: nodeNr as index, array according to |nodeNr of pos-sibling|nodeNr of neg-sibling|; -1 as 
    //                         nodeNr of a sibling indicates no sibling at this position. An example:
    //                                    0             [0]->| 1| 2|
    //                                  /   \      =>   [1]->|-1|-1|
    //                                 1     2          [2]->| 4| 5|
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
    SC_ClassifierTree::SC_ClassifierTreeNode* constructTree(int nodeCount, int **treeConstructionTable, int **posLabels, int *posLabelsCount, int **negLabels, int *negLabelsCount, int *posResult, int *negResult, unsigned int classifierType);

    //====================================================================================================================
    //	this method corrects an error that may have been occured during training: if in the tree-construction table 
    //  positive- and negative nodes where switched, this can be corrected with this method for an already trained, but
    //  wrongly connected classifier-tree
    //====================================================================================================================
    void switchLeftRightNodes(SC_ClassifierTree::SC_ClassifierTreeNode *pNode);

    //====================================================================================================================
    //	returns true if inconsistencies seem to exist in the connection of the nodes according to the results of parents
    //  and theire siblings
		//  TODO: assumes result-types to be bitflags in 32bit variables
    //====================================================================================================================
    bool resultContradiction(SC_ClassifierTree::SC_ClassifierTreeNode *pNode);

  public :

    SC_ClassifierTree(SC_TweakableParameters* pTweak, bool doScaling, bool verbose = true, bool traverseCompleteTree = false, int nodeCount = 0, int **treeConstructionTable = NULL, int **posLabels = NULL, int *posLabelsCount = NULL, int **negLabels = NULL, int *negLabelsCount = NULL, int *posResult = NULL, int *negResult = NULL, unsigned int classifierType = sclib::ctSVM);
    virtual ~SC_ClassifierTree();

    //====================================================================================================================
    //	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
    //====================================================================================================================
    virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative);

    //====================================================================================================================
    //	train a classifier for distinguishing between several classes
    //  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
    //  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
    //  respective row of pData.
    //====================================================================================================================
    virtual int trainMultiClass(SV_Data *pData, int *classes);

		//====================================================================================================================
		//	classifiy previously unseen test-data; returned is an array of classlabels, each entry corresponding to the 
		//  respective row in pData; if available, the probabilities for each class-decision are given in the pProbabilities 
		//  parameter, but here, the array-indexing has a different meaning than in the other classifier classes:
		//  in col. nodeNr*2 is the label/labels assigned by this node, in nodeNr*2+1 is the correspondig probability
		//
		//  ATTENTION: only a previously trained (or loaded) classifier can be used for classification!
		//====================================================================================================================
    virtual int* classify(SV_Data *pData, SV_Data* &pProbabilities);

    //====================================================================================================================
    //	save and load a trained classifier to/from a file
    //====================================================================================================================
    virtual int saveClassifier(const char *fileName);
    virtual int loadClassifier(const char *fileName);

    //====================================================================================================================
    //	to convert between given labels and indices into the probability-parameter of the classifiy()-method
    //====================================================================================================================
		virtual long int label2idx(long int label);
		virtual long int idx2label(long int idx);

    //====================================================================================================================
    //	returns the type of the used sub-classifier (type of the nodes in the tree)
    //====================================================================================================================
		int getClassifierType(void) {return this->classifierType;}
};

#endif
