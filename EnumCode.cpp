#include "EnumCode_v2.h"
#include <cstdlib>
#include <vector>
#include <gmp.h>

using namespace std;

mpz_t binomCoeff[MAXBLOCKLEN][MAXBLOCKLEN+1];

//Constructors-----------------------------------------------------------------------------------

void Helper_Constructor( )
{
    unsigned int i,j;
    for(i=0; i<MAXBLOCKLEN; i++)
    {
        mpz_init(binomCoeff[i][0]);
        mpz_set_ui(binomCoeff[i][0],1);
        for(j=1; j<=i; j++)
        {
            mpz_init(binomCoeff[i][j]);
            mpz_add(binomCoeff[i][j],binomCoeff[i-1][j-1],binomCoeff[i-1][j]);
        }
        mpz_init(binomCoeff[i][j]);
    }
}

EnumCode::EnumCode( )
{
	Helper_Constructor();
	alphabet = new unsigned char[2];
}

EnumCode::EnumCode( unsigned alphabetSize )
{
	Helper_Constructor();
    this->alphabetSize = alphabetSize;    
	alpVector.resize(alphabetSize);
	alphabet = new unsigned char[alphabetSize];
}

EnumCode::EnumCode( vector<unsigned int> v ) 
{
	Helper_Constructor();
    this->alphabetSize = v.size();
    this->alpVector = v;
	sequence.resize(GetInnerSum());	
	alphabet = new unsigned char[alphabetSize];
}

EnumCode::EnumCode( vector<unsigned int> v, vector<unsigned int> s )
{
	Helper_Constructor();
    this->alphabetSize = v.size();
	this->alpVector = v;	
	this->sequence = s;
	alphabet = new unsigned char[alphabetSize];
}

EnumCode::~EnumCode( )
{
	alpVector.clear();
	sequence.clear();
	delete[] alphabet;
}


//Set function members------------------------------------------------------------------------

void EnumCode::SetAlphaSize( unsigned alphabetSize )
{
	this->alphabetSize = alphabetSize;
	alpVector.resize(alphabetSize);
	for (unsigned int i=0; i<alphabetSize; i++)
		alpVector[i] = 0;
}

void EnumCode::SetVector( vector<unsigned int> v )
{
	SetAlphaSize(v.size());
	this->alpVector = v;
	sequence.resize(GetInnerSum());	
}

void EnumCode::SetSequence( vector<unsigned int> sequence )
{
	unsigned int i, max=0;
	this->sequence = sequence;
	
	for (i=0; i<sequence.size(); i++)
		if (max < sequence[i])
			max = sequence[i];

	SetAlphaSize(max);
    for (i=0; i<sequence.size(); i++)
		alpVector[sequence[i]-1]++;
}

void EnumCode::SetAlphabet( unsigned char* alpha )
{
	delete[] alphabet;
	alphabet = new unsigned char[alphabetSize];
	for (unsigned int i=0; i<alphabetSize; i++)
		alphabet[i] = alpha[i];
}


//Get function members------------------------------------------------------------------------

unsigned int EnumCode::GetAlphaSize( )
{
	return this->alphabetSize;	
}

vector<unsigned int> EnumCode::GetVector( )
{
	return this->alpVector;
}

vector<unsigned int> EnumCode::GetSequence( )
{
	return this->sequence;
}

unsigned int EnumCode::GetInnerSum( )
{
	unsigned int i, innerSum = 0;
    for (i = 0; i<GetAlphaSize(); i++)
		innerSum += alpVector[i];	
	return innerSum;
}

unsigned char* EnumCode::GetAlphabet( )
{
	return this->alphabet;
}


//--------------------------------------------------------------------------------------------
//functions 
void EnumCode::Number_to_BinarySeq( unsigned int number )
{
	unsigned int i=0, temp=number;

	while (temp > 1)
	{
		temp /= 2;	
		i++;	
	}	

	sequence.resize(i);
	i--;
	while (number > 1)
	{
		sequence[i--] = number%2 + 1;
		number /= 2;
	}	
	SetSequence(sequence);
}


void EnumCode::Number_to_BinarySeq( mpz_t number )
{
	unsigned int i=0;
    mpz_t temp;
	mpz_init_set(temp, number);

	while (mpz_cmp_ui(temp,1) > 0)
	{
		mpz_divexact_ui(temp,temp,2);	
		i++;	
	}	

	sequence.resize(i);
	i--;
	
	mpz_set(temp, number);
	while (mpz_cmp_ui(temp,1) > 0)
	{
		if (mpz_odd_p(temp))
			sequence[i--] = 2;
		else
			sequence[i--] = 1;

		mpz_divexact_ui(temp,temp,2);
	}	

	SetSequence(sequence);
    mpz_clear(temp);	
}


//---------------------------------------

unsigned EnumCode::BinarySeq_to_Number( )
{	
	int i;	
	unsigned result;    
	mpz_t total, op1, op2;
	mpz_init_set_ui(total, 0);	
	mpz_init_set_ui(op2, 1);
	mpz_init(op1);
	
	if (sequence.empty())
		return 0;
	
	for (i=sequence.size()-1; i>=0; i--)
	{
		mpz_set_ui(op1,sequence[i]-1);
		mpz_addmul(total,op1,op2);					
		mpz_mul_ui(op2,op2,2);
	}

	mpz_add(total,total,op2);
	result = mpz_get_ui(total);
	
    mpz_clear(op1);
    mpz_clear(op2);
    mpz_clear(total);

	return result;
}


unsigned EnumCode::BinarySeq_to_Number( mpz_t total )
{
    mpz_t op1, op2;
	mpz_set_ui(total, 0);	
	mpz_init_set_ui(op2, 1);
	mpz_init(op1);
		
	if (sequence.empty())
		return 0;

	for (int i=sequence.size()-1; i>=0; i--)
	{
		mpz_set_ui(op1,sequence[i]-1);
		mpz_addmul(total,op1,op2);					
		mpz_mul_ui(op2,op2,2);
	}

	mpz_add(total,total,op2);
    mpz_clear(op1);
    mpz_clear(op2);

	return 1;
}


//-----------------------------------------

int EnumCode::Index_to_Vector_ui( unsigned block_length, unsigned counter, unsigned index )
{
    mpz_t temp, temp2;
	mpz_init_set_ui(temp, index);	
	mpz_init(temp2);

	SetAlphaSize(counter);
    Count_ParVectors(block_length, alphabetSize, temp2);
    if (mpz_cmp(temp,temp2) >= 0)
		return -1; 

	Helper_Index_to_Vector(block_length, temp);
    mpz_clear(temp);
    mpz_clear(temp2);	

	return 0;
}


int EnumCode::Index_to_Vector_ui( unsigned index )
{
	unsigned int block_length = GetInnerSum();
    mpz_t temp, temp2;
	mpz_init_set_ui(temp, index);	
	mpz_init(temp2);

    Count_ParVectors(block_length, alphabetSize, temp2);
    if (mpz_cmp(temp,temp2) >= 0)
		return -1; 

	Helper_Index_to_Vector(block_length, temp);
    mpz_clear(temp);
    mpz_clear(temp2);	

	return 0;
}

int EnumCode::Index_to_Vector( unsigned block_length, unsigned counter, mpz_t index )
{
	mpz_t temp;
	mpz_init(temp);

	SetAlphaSize(counter);
    Count_ParVectors(block_length, alphabetSize, temp);
    if (mpz_cmp(index,temp) >= 0)
		return -1; 

	Helper_Index_to_Vector(block_length, index);
    mpz_clear(temp);

	return 0;
}


int EnumCode::Index_to_Vector( mpz_t index )
{
	unsigned int block_length = GetInnerSum();
	mpz_t temp;
	mpz_init(temp);

    Count_ParVectors(block_length, alphabetSize, temp);
    if (mpz_cmp(index,temp) >= 0)
		return -1; 

	Helper_Index_to_Vector(block_length, index);
    mpz_clear(temp);

	return 0;
}


//-----------------------------------------

unsigned EnumCode::Vector_to_Index( )
{
    unsigned int block_length = GetInnerSum();
    unsigned int dimension, i;
    mpz_t index, temp;
    mpz_init(temp);
	mpz_init_set_ui(index,0);

	if (alpVector.empty())
		return 0;

    for (dimension = 0; dimension<(alphabetSize-1); dimension++)
    {
        for (i=0; i<alpVector[dimension]; i++)
        {
            Count_ParVectors(block_length-i, alphabetSize-dimension-1, temp);
            mpz_add(index, index, temp);
        }
        block_length -= alpVector[dimension];
    }
	i = mpz_get_ui(index);
	
    mpz_clear(index);
    mpz_clear(temp);

	return i;
}


unsigned EnumCode::Vector_to_Index( mpz_t index )
{
    unsigned int block_length = GetInnerSum();
    unsigned int dimension, i;
    mpz_t temp;
    mpz_init(temp);
	mpz_set_ui(index,0);

	if (alpVector.empty())
		return 0;

    for (dimension = 0; dimension<(alphabetSize-1); dimension++)
    {
        for (i=0; i<alpVector[dimension]; i++)
        {
            Count_ParVectors(block_length-i, alphabetSize-dimension-1, temp);
            mpz_add(index, index, temp);
        }
        block_length -= alpVector[dimension];
    }
    mpz_clear(temp);

	return 1;
}


//-----------------------------------------

unsigned EnumCode::PermIndex_to_Seq( mpz_t permID )
{
	mpz_t temp;
	mpz_init(temp);	

	if (alpVector.empty())
		return 0;

	NumberofPermutations(temp);
	if (mpz_cmp(permID, temp) >= 0)
		return 0;

	Helper_PermIndex_to_Seq(permID);	
	
	return 1;	
}


//-----------------------------------------

unsigned EnumCode::Seq_to_PermIndex( )
{
	unsigned int i, vectorSize;
    mpz_t permID;
    mpz_init_set_ui(permID, 0);
 
	if (sequence.empty())
		return 0;

    for (i=0; i<sequence.size(); i++)
		if (vectorSize < sequence[i])
			vectorSize = sequence[i];
	Helper_Seq_to_PermIndex	(permID, vectorSize); 
	i = mpz_get_ui(permID);

    mpz_clear(permID);
	return i;
}


unsigned EnumCode::Seq_to_PermIndex( mpz_t permID )
{
	unsigned int i, vectorSize;

	if (sequence.empty())
		return 0;

    for (i=0; i<sequence.size(); i++)
		if (vectorSize < sequence[i])
			vectorSize = sequence[i];
	Helper_Seq_to_PermIndex	(permID, vectorSize); 

	return 1;
}


//-----------------------------------------

/*
unsigned EnumCode::Count_ParVectors( unsigned int sum, unsigned int counter )
{
	unsigned temp;
    mpz_t op1, op2, op3;
    mpz_init_set_ui(op1, sum-1);
	mpz_add_ui(op1, op1, counter);		
	mpz_init_set(op3, op1);
    mpz_init_set_ui(op2, 2);
	
	temp = min(sum, counter)-1;
	for (counter=temp; counter>1; counter--)
	{
		mpz_sub_ui(op3, op3, 1);
		mpz_mul(op1, op1, op3);
		while (mpz_divisible_p(op1, op2) && temp>1)
		{
			mpz_divexact(op1, op1, op2);
			mpz_add_ui(op2, op2, 1);
			temp--;
		}
	}
	temp = mpz_get_ui(op1);
    mpz_clear(op1);
    mpz_clear(op2);	
	mpz_clear(op3);

	return temp;
}
*/

unsigned EnumCode::Count_ParVectors( unsigned int sum, unsigned int counter )
{
	unsigned int i;	
	mpz_t value;
	mpz_init(value);	
	
    mpz_t op1, op2;
    mpz_init(op1);
    mpz_init(op2);
    mpz_set_ui(value,0); 
	if (sum<1)
		return 0;

    for(i=0; i<counter; i++)
    {
        mpz_set(op1,binomCoeff[sum-1][counter-(i+1)]);
        mpz_set(op2,binomCoeff[counter][i]);
        mpz_addmul(value,op1,op2);
    }
	i = mpz_get_ui(value);

    mpz_clear(op1);
    mpz_clear(op2);
    mpz_clear(value);	

	return i;
}

/*
unsigned EnumCode::Count_ParVectors( unsigned int sum, unsigned int counter, mpz_t op1 )
{
	unsigned temp;
    mpz_t op2, op3;
    mpz_set_ui(op1, sum-1);
	mpz_add_ui(op1, op1, counter);		
	mpz_init_set(op3, op1);
    mpz_init_set_ui(op2, 2);
	
	temp = min(sum, counter)-1;
	for (counter=temp; counter>1; counter--)
	{
		mpz_sub_ui(op3, op3, 1);
		mpz_mul(op1, op1, op3);
		while (mpz_divisible_p(op1, op2) && temp>1)
		{
			mpz_divexact(op1, op1, op2);
			mpz_add_ui(op2, op2, 1);
			temp--;
		}
	}
    mpz_clear(op2);	
	mpz_clear(op3);

	return 1;
}
*/

unsigned EnumCode::Count_ParVectors( unsigned int sum, unsigned int counter, mpz_t value )
{
    mpz_t op1, op2;
    mpz_init(op1);
    mpz_init(op2);
    mpz_set_ui(value,0);
    if (sum<1)
		return 0;

    for(unsigned int i=0; i<counter; i++)
    {
        mpz_set(op1,binomCoeff[sum-1][counter-(i+1)]);
        mpz_set(op2,binomCoeff[counter][i]);
        mpz_addmul(value,op1,op2);
    }
    mpz_clear(op1);
    mpz_clear(op2);

	return 1;
}


//--------------------------------------------------------------------------------------------

unsigned EnumCode::Number_to_BinaryIndex( unsigned int num )
{	
	unsigned int i;
	mpz_t index;
	mpz_init(index);		
	Number_to_BinarySeq(num);

	if (sequence.empty())
		return 0;

	Seq_to_PermIndex(index);
 	i = mpz_get_ui(index);
    mpz_clear(index);	
	return i;
}


unsigned EnumCode::Number_to_BinaryIndex( unsigned int num, mpz_t index )
{
	Number_to_BinarySeq(num);

	if (sequence.empty())
		return 0;		

	Seq_to_PermIndex(index);

	return 1;
}


//--------------------------------------------------------------------------------------------

unsigned EnumCode::BinaryIndex_to_Number( unsigned int index, vector<unsigned int> vec )
{
	unsigned int i;
	alpVector.resize(2);
	alpVector[0] = vec[1];
	alpVector[1] = vec[0]-vec[1];

	mpz_t temp;
	mpz_init_set_ui(temp, index);
	PermIndex_to_Seq(temp);
	BinarySeq_to_Number(temp);

 	i = mpz_get_ui(temp);
    mpz_clear(temp);	
	return i;
}


void EnumCode::BinaryIndex_to_Number( unsigned index, vector<unsigned int> vec, mpz_t value )
{
	alpVector.resize(2);
	alpVector[0] = vec[1];
	alpVector[1] = vec[0]-vec[1];

	mpz_set_ui(value, index);
	PermIndex_to_Seq(value);
	BinarySeq_to_Number(value);
}


//--------------------------------------------------------------------------------------------
//Private functions
void EnumCode::Helper_Index_to_Vector( unsigned int block_length, mpz_t index )
{
    unsigned int dimension = 0;
    mpz_t temp, tempIndex;
	mpz_init(temp);
    mpz_init_set(tempIndex, index);

    for (dimension = 0; dimension<(alphabetSize-1); dimension++)
    {
        alpVector[dimension] = 0;
        while (block_length > 0)
        {
            Count_ParVectors(block_length, alphabetSize-dimension-1, temp);
            if (mpz_cmp(tempIndex,temp) >= 0)
            {
                mpz_sub(tempIndex, tempIndex, temp);
                alpVector[dimension]++;
                block_length--;
            }
            else break;
        }
    }
	    
	alpVector[dimension] = block_length;
	mpz_clear(temp);
    mpz_clear(tempIndex);	
}


void EnumCode::Helper_PermIndex_to_Seq( mpz_t permID )
{
    unsigned int i, j, block_length = GetInnerSum();
   	vector<unsigned int> Vector = alpVector;
    mpz_t pid, t, temp;
    mpz_init_set(pid, permID);
    mpz_init(t);
    mpz_init(temp);

	sequence.resize(block_length);

    for (i=0; i<block_length; i++)
    {
        sequence[i] = 1;
        mpz_set_ui(temp,0);
        for (j=0; j<alphabetSize; j++)
            if (alpVector[j] > 0)
            {
                alpVector[j]--;
                NumberofPermutations(t);
                mpz_add(temp,temp,t);

                if (mpz_cmp(pid,temp) >= 0)
                    alpVector[j]++;
                else
                {
                    sequence[i] = j+1;
                    mpz_sub(temp,temp,t);
                    mpz_sub(pid,pid,temp);
                    break;
                }
            }
    }

	alpVector = Vector;
}


void EnumCode::Helper_Seq_to_PermIndex( mpz_t permID, unsigned AlphaSize )
{
    unsigned int i, j;
	mpz_set_ui(permID, 0);
    mpz_t t;
    mpz_init(t);	
	SetAlphaSize(AlphaSize);
	
    for (i=0; i<sequence.size(); i++)
		alpVector[sequence[i]-1]++;

    for (i=0; i<sequence.size(); i++)
    {
        for (j=0; j<sequence[i]-1; j++)
        {
            if (alpVector[j]>0)
            {
                alpVector[j]--;
                NumberofPermutations(t);
                mpz_add(permID, permID, t);
                alpVector[j]++;
            }
        }
        alpVector[sequence[i]-1]--;
    }
}


void EnumCode::NumberofPermutations( mpz_t value )
{
    mpz_t denominator, nominator;
    mpz_init_set_ui(denominator,1);
    mpz_init_set_ui(nominator  ,1);
    unsigned int i,j,total=0;

    for (i=1; i<=GetInnerSum(); i++)
		mpz_mul_ui(nominator,nominator,i);

    for (i=0; i<GetAlphaSize(); i++)
        for (j=1; j<=alpVector[i]; j++)
            mpz_mul_ui(denominator, denominator, j);

    mpz_divexact(value, nominator, denominator);
    mpz_clear(nominator);
    mpz_clear(denominator);
}


//--------------------------------------------------------------------------------------------

void message( )
{
        cout << " 1. Count Distinct Parikh Vectors     : Find how many distinct ways can K positive integers sum up to S " << endl;
        cout << " 2. List  All Parikh Vectors          : Find how many distinct ways can K positive integers sum up to S and list all of them" << endl;
        cout << " 3. Parikh vetor to Index             : Find the enumerative index of given frequencies of alphabet characters" << endl;
        cout << " 4. Index to Parikh vector            : Given the enumerative index and the length of the sequence (inner sum of Parikh vector), find the frequencies of alphabet characters" << endl;
        cout << " 5. Number of permutations            : Given the frequencies of alphabet characters, find the number of distinct permutations" << endl;
        cout << " 6. List all permutations             : Given the frequencies of alphabet characters, list all distinct permutations" << endl;
        cout << " 7. PermIndex to Sequence             : Given the frequencies of alphabet characters, and the permutation index, retrieve the corresponding permutation" << endl;
        cout << " 8. Sequence to PermIndex             : Given a permutation of alphabet characters, find the corresponding permutation index" << endl;
        cout << " 9. Binary Sequence to Index          : Given a number, find binary sequence and its permutation index" << endl;
	    cout << "10. Index to Binary Sequence          : Given a frequencies vector and a index, find the corresponding number and its binary sequence" << endl;	
        cout << "11. Exit"  << endl;		
        cout << "ENTER COMMAND ID: \t" ;
}


void EnumCode::test()
{
    unsigned int command, K, blockLength, a, i;
    char symbol, tmpstr[1000];
    unsigned char pattern[1024]= {0};
	unsigned char temp[1024];
	vector <unsigned int> vec;
    mpz_t index,tempmpz;
    mpz_init(index);
    mpz_init(tempmpz);

    while (1)
    {
		message();
        scanf("%d", &command);
        switch (command)
        {

        case 1:
            cout << "Enter K (number of counters):\t";
            scanf("%d",&K);
            cout << "Enter Sum (sum of counters): \t";
            scanf("%d", &blockLength);
            Count_ParVectors(blockLength,K,tempmpz);
            memset(tmpstr,0,1000);
            mpz_get_str(tmpstr,10,tempmpz);
            cout << "Number of distinct " << K << " dimensional parikh vectors having an inner sum of " << blockLength << " is:\t" << tmpstr<< endl<<endl;
            break;

        case 2:
            cout << "Enter K (number of counters):\t";
			scanf("%d",&K);
            cout << "Enter Sum (sum of counters): \t";
            scanf("%d", &blockLength);
            SetAlphaSize(K);
            Count_ParVectors(blockLength,K,tempmpz);
            memset(tmpstr,0,1000);
            mpz_get_str(tmpstr,10,tempmpz);
            cout << "Number of distinct " << K << " dimensional parikh vectors having an inner sum of " << blockLength <<" is " << tmpstr<< endl;
            cout << "Do you want to list all such vectors (y/n)?\t";
            scanf("%c",&symbol);
            scanf("%c",&symbol);
            if (symbol == 'y')
            {
                mpz_set_ui(index,0);
                while(mpz_cmp(index,tempmpz)<0)
                {
                    Index_to_Vector(blockLength,K,index);
                    cout << mpz_get_str(tmpstr,10,index);
                    for(a=0; a<K; a++) cout << '\t' << GetVector()[a];
                    cout << endl;
                    mpz_add_ui(index,index,1);
                }
            }
            cout << endl << endl;
            break;

        case 3:
            cout << "What is the dimension of the Parikh vectors (number of counters):\t";
            scanf("%d",&K);
            SetAlphaSize(K);
            cout << "Enter the values for each dimension one by one (e.g., 0 1 2 1):  \t";

			vec.resize(K);
            for(a=0,i=0; i<K; i++)
            {
                scanf("%d",&vec[i]);
                a+=vec[i];
            }
			SetVector (vec);
            Vector_to_Index(tempmpz);
            cout << "This vector corresponds to index " << mpz_get_str(tmpstr,10,tempmpz) << " in the list of all possible vectors  summing up to " << a  << endl << endl;
            break;

        case 4:
            cout << "What is the dimension of the Parikh vectors (number of counters):\t";
            scanf("%d",&K);
            SetAlphaSize(K);
            cout << "Enter block length (inner sum of the Parikh vector):\t";
            scanf("%d", &blockLength);
            cout << "Enter index        ";
            scanf("%s", tmpstr);
            mpz_set_str(tempmpz,tmpstr,10);
            if (-1 != Index_to_Vector(blockLength,K,tempmpz))
            {
                cout << "Corresponding vector is " << '\t';
                for(i=0; i<K; i++) cout << GetVector()[i] << "  ";
                cout << endl << endl;
            }
            else
            {
                cout << "Index is invalid" << endl << endl;
            }
            break;

        case 5:
            cout << "What is the dimension of the Parikh vectors (number of counters):\t";
            scanf("%d",&K);
            SetAlphaSize(K);
            cout << "Enter the values for each dimension one by one (e.g., 0 1 2 1):  \t";
			vec.resize(K);
            for(a=0,i=0; i<K; i++)
            {
                scanf("%d",&vec[i]);
                a+=vec[i];
            }
			SetVector (vec);
            NumberofPermutations(tempmpz);
            cout << "Number of possible distinct permutations is \t" << mpz_get_str(tmpstr,10,tempmpz) <<endl << endl;
            break;

        case 6:
            cout << "What is the dimension of the Parikh vectors (number of counters):\t";
            scanf("%d",&K);
            SetAlphaSize(K);
            cout << "Enter the alphabet without spaces (e,g,. abcd): \t";
            scanf("%s",temp);
			SetAlphabet(temp);
            cout << "Enter the values for each dimension one by one (e.g., 0 1 2 1): \t";
			vec.resize(K);
            for(a=0,i=0; i<K; i++)
            {
                scanf("%d",&vec[i]);
                a+=vec[i];
            }
			SetVector (vec);
            NumberofPermutations(tempmpz);
            cout << "Number of possible distinct permutations is\t" << mpz_get_str(tmpstr,10,tempmpz) << endl;
            cout << "Do you really want to list all these permutations (y/n)? \t";
            scanf("%c",&symbol);
            scanf("%c",&symbol);
            if (symbol == 'y')
            {
                mpz_set_ui(index,0);
                while(mpz_cmp(index,tempmpz)<0)
                {
                    memset(pattern,0,1000);
                    PermIndex_to_Seq(index);
                    cout << mpz_get_str(tmpstr,10,index) << '\t' ;
                    for(i=0; i<a; i++) cout << GetAlphabet()[(GetSequence()[i]-1)] << "  ";
                    cout << endl;
                    mpz_add_ui(index,index,1);
                }
            }
            cout <<endl<<endl;
            break;

        case 7:
            cout << "What is the dimension of the Parikh vectors (number of counters):\t";
            scanf("%d",&K);
            SetAlphaSize(K);
            cout << "Enter the alphabet without spaces (e,g,. abcd):  \t";
            scanf("%s",temp);
			SetAlphabet(temp);
            cout << "Enter the values for each dimension one by one (e.g., 0 1 2 1): \t";
			vec.resize(K);
            for(a=0,i=0; i<K; i++)
            {
                scanf("%d",&vec[i]);
                a+=vec[i];
            }
			SetVector (vec);
            cout << "Enter the index:\t";
            scanf("%s",tmpstr);
            mpz_set_str(tempmpz,tmpstr,10);
            memset(pattern,0,1000);
            if (-1!=PermIndex_to_Seq(tempmpz))
            {
                cout << "The corresponding permutation of the input index is \t" ;
                    for(i=0; i<a; i++) cout << GetAlphabet()[(GetSequence()[i]-1)] << "  ";
                cout << endl<<endl;
            }
            else
            {
                cout<< "Invalid index " << endl << endl;
            }
            break;

        case 8:
            cout << "Enter the alphabet size:\t";
            scanf("%d",&K);
            SetAlphaSize(K);

            cout << "Length of the sequence (inner sum of Parikh vector):\t";
            scanf("%d",&a);

            cout << "Enter a permutation via sequence of integers (e.g., 2 1 1 4 2 3 ):\t" ;
			vec.resize(a);
            for(i=0; i<a; i++)
                scanf("%d",&vec[i]);
			SetSequence (vec);

            Seq_to_PermIndex(tempmpz);
            cout << "The input permutation corresponds to index " << mpz_get_str(tmpstr,10,tempmpz) << " in the list of all permutations with respective symbol occurrences. " << endl<<endl;
            break;

		case 9:
			cout << "Enter the number:\t";		
			scanf("%d",&K);
			Number_to_BinarySeq(K);
			
			cout << "The corresponding binary sequence of the input number is (missing the first one):\t" ;
            for(i=0; i<GetSequence().size(); i++) cout << GetSequence()[i]-1 << "  ";
                cout << endl;
			
			Seq_to_PermIndex(tempmpz);	
			cout << "The input permutation corresponds to index " << mpz_get_str(tmpstr,10,tempmpz) << " in the list of all permutations with respective symbol occurrences. " << endl<<endl;
			break;

        case 10:
			cout << "Enter the index:\t";		
            scanf("%s",tmpstr);
            mpz_set_str(index,tmpstr,10);		
            SetAlphaSize(2);

			vec.resize(2);
			cout << "Enter the length of the binary sequence:\t";
			scanf("%d",&vec[1]);

            cout << "Enter the length of the binary sequence and the number of 0's:\t";
			scanf("%d",&vec[0]);

			vec[1] -= vec[0];
			SetVector (vec);

			PermIndex_to_Seq(index);
			cout << "The corresponding binary sequence of the input number is (missing the first one)" << '\t' ;
            for(i=0; i<GetSequence().size(); i++) cout << GetSequence()[i]-1 << "  ";
                cout << endl;

			BinarySeq_to_Number(tempmpz);
			cout << "The number is " << mpz_get_str(tmpstr,10,tempmpz) << endl<<endl;
			break;	
			
		case 11:
            return;
        }
    }

}
