#ifndef ENUMCODE_H
#define ENUMCODE_H
#include <iostream>
#include <memory.h>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <gmp.h>

#define MAXBLOCKLEN 1024

using namespace std;

class EnumCode
{	
    public:
		//Constructors
		EnumCode	();
		EnumCode	(unsigned alphabetSize);
        EnumCode 	(vector<unsigned int> alpVector);
		EnumCode	(vector<unsigned int> alpVector, vector<unsigned int> sequence);        
		~EnumCode	();
     
		//Set function members   
		void SetAlphaSize 	(unsigned alphabetSize); 
		void SetVector 		(vector<unsigned int> alpVector); 
		void SetSequence	(vector<unsigned int> sequence);
		void SetAlphabet 	(unsigned char* alphabet); 				
		
		//Get function members
		unsigned int			GetAlphaSize();
		vector<unsigned int> 	GetVector 	();
		vector<unsigned int> 	GetSequence	();		 		
		unsigned int			GetInnerSum	();
		unsigned char*			GetAlphabet ();
        
		//functions    
		void 	 Number_to_BinarySeq	(unsigned number);  	
		void 	 Number_to_BinarySeq	(mpz_t number); 

		unsigned BinarySeq_to_Number	(); 
		unsigned BinarySeq_to_Number	(mpz_t value); 				
		
		int		 Index_to_Vector_ui 	(unsigned sum, unsigned counter, unsigned index); 
		int		 Index_to_Vector_ui 	(unsigned index);
		int		 Index_to_Vector 		(unsigned sum, unsigned counter, mpz_t index);
		int		 Index_to_Vector 		(mpz_t index);

        unsigned Vector_to_Index 		(); 
        unsigned Vector_to_Index 		(mpz_t index);
        
		unsigned PermIndex_to_Seq	 	(mpz_t permID);

        unsigned Seq_to_PermIndex		();
        unsigned Seq_to_PermIndex		(mpz_t permID); 

		unsigned Count_ParVectors 		(unsigned sum, unsigned counter);
		unsigned Count_ParVectors 		(unsigned sum, unsigned counter, mpz_t value); 
		
		unsigned Number_to_BinaryIndex	(unsigned num); 
		unsigned Number_to_BinaryIndex	(unsigned num, mpz_t value);
	
		unsigned BinaryIndex_to_Number	(unsigned index, vector<unsigned int> vec);
		void 	 BinaryIndex_to_Number	(unsigned index, vector<unsigned int> vec, mpz_t value);
		
		void test();

    private:
        unsigned int  				alphabetSize;
        vector<unsigned int> 		alpVector;
		vector<unsigned int> 		sequence;
		unsigned char*				alphabet;

		void Helper_Index_to_Vector 	(unsigned block_length, mpz_t index);	
		void Helper_PermIndex_to_Seq 	(mpz_t permID);
        void Helper_Seq_to_PermIndex	(mpz_t permID, unsigned vectorSize); 
		void NumberofPermutations 		(mpz_t value);	
};

#endif // ENUMDNA_H
