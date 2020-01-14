/*  AMAN
 *  KHURANA
 *  AMANKHUR
 */

#ifndef A0_HPP
#define A0_HPP
#include<iostream>
#include<vector>
#include<omp.h>
#include<math.h>
#include <typeinfo>


using namespace std;

template <typename T,typename Op>
// this function performs reduction of a given input vector
// and an operator
// returns the reduced array
vector<T> reduce(vector<T> a, Op op, T identity_element){


// size of input to reduce
int n = a.size();
// extracting closest exponent of 2
// used for determing padding length
int max_pow = floor(log2(n));

// calculating padding length
int pad_length = pow(2,max_pow+1) - n;

//int identity_element = 0;

//creating vector of padding length filled with the identity element 
vector<T> padding(pad_length,identity_element);

//doing padding
a.insert(a.end(),padding.begin(),padding.end());

// new n 
int n_effective = n+ pad_length;



// iterating through levels of our conceptual binary tree
// and performing reduction
// taking right and left child of the root, applying the
// operator on them and updating the value of the right child
// with the calculated value

for (int d=0; d < log2(n_effective); d++){
        int inc = pow(2,d+1); 
        #pragma omp parallel for
        for (int i=0; i<n_effective; i= i + inc){
            a[i+pow(2,d+1)- 1] = op(a[i+pow(2,d) - 1],a[i+pow(2,d+1)-1]);
        }
    }

return a;

}

template <typename T, typename Op>
// this function takes as input the reduced array from the "reduce"
// function and returns prescan of the original array which is the input
// of "reduce".
// n_original is required to slice the array incase of extra padding
// to take the relevant part of the array as the output
//
vector<T> downSweep(vector<T> a, int n_original, int a_last, Op op, T identity_element){
    
    int n = a.size();
    //int identity_element = 0;

    // assigning the last element as identity
    a[n-1] = identity_element;

   //  Here we iterate down in our conceptual binary tree
   //  to calculate prescan
    for (int d= log2(n)-1; d>=0; d--){
        
        int inc = pow(2,d+1); 
#pragma omp parallel for
        for(int i=0; i<n; i = i + inc){
          
            // storing the value of the left child 
            int temp = a[i+pow(2,d)- 1]; 
            
            // assigning value of right child to the left child 
            a[i+pow(2,d)- 1] = a[i+pow(2,d+1)- 1]; 
          
            // assigning the value of the operator applied to the value of 
            // left child stored in temp and the right child 
            a[i+pow(2,d+1)- 1] = op(temp,a[i+pow(2,d+1)- 1]);

        }

    } 
    vector<T> prescan = vector<T>(a.begin(),a.begin()+n_original);

return prescan;
}

template <typename T, typename Op>
void omp_scan(int n, const T* in, T* out, Op op) {

    // buffering original n to get correct result in case of paddings
    int n_original = n;

    // a is our input vector
    vector<T> a(in,in+n);

    // vector to store intermediate output calculations,
    // we use this to get the final result
    vector<T> output(n);

    // storing the last element to convert prescan to scan
    int last_element = a[n-1];

    // vector to store processor sums (bad naming!)
    vector<T> sum;
    
    // p_max will store the number of available processors
    int p_max = 0;
    int n_mod_p = 0;

    // stores the identity element
    string op_name = typeid(op).name();

    // extracting identity element from operator type   

    T identity_element = 0 ;
    string plus = "plus";
    string mul = "multi";
 
    if(op_name.find(plus) != std::string::npos){
    identity_element = 0;
}
    if(op_name.find(mul) != std::string::npos){
    identity_element = 1;
}


#pragma omp parallel
    {
        #pragma omp single 
        {
        p_max = omp_get_num_threads(); // initializing p_max
        sum.resize(p_max);  // initalizing processor sums vector
        n_mod_p = n % p_max;
	
        // to handle the values in case when n/p is not an integer
        if(n_mod_p !=0){
                vector<T> padding(p_max - n_mod_p,identity_element);
                a.insert(a.end(),padding.begin(),padding.end());
		n = a.size();
                last_element = a[n-1];
		output.resize(n);
		}

	}    
        #pragma omp for
        // for each processor we assign a chunk of the input
        // and calculate the total sum for each processor
        // and store the resulting sum for each processor in
        // the vector "sum" 
            for(int i=0;i<p_max; i++){
                int subarray_length =  floor(n/p_max);
                sum[i] = a[subarray_length*i];
                T i_sum = sum[i];
                #pragma omp parallel for reduction(+:i_sum)
                for(int j=1; j <subarray_length; j++)
                    i_sum = op(i_sum, a[subarray_length*i + j]);
                
                sum[i] = i_sum;
        }
    }

 
    // we perform reduction and downsweep on the processor
    // sums, we will use these sums to retrieve the final prescan
    // of the original input from results for our individual chunks    
    vector<T> result = reduce(sum,op,identity_element);
    vector<T> sum_prescan = downSweep(result,sum.size(),sum[sum.size()-1],op, identity_element);
    

   #pragma omp parallel
   {
    #pragma omp single
    {
    p_max = omp_get_num_threads();
    }   
       
        // for each processor we assign a chunk of the input
        // then we perform reduction and downsweep to each of the input
        // chunk.
        // finally we use the individual results of the chunks along 
        // with the processor sums array(sum) to get the prescan of the 
        // original input  
     
    #pragma omp for
        for(int i=0;i<p_max; i++){

	    int subarray_length =  floor(n/p_max);
            vector<T> subarray(subarray_length);
            subarray[0] = a[subarray_length*i];
            
            // extracting subarray for each of the
            // i'th processor
            for(int j=1; j <subarray_length; j++)
                subarray[j] = a[subarray_length*i + j];

            
            
            vector<T> subarray_scan = downSweep(reduce(subarray,op,identity_element),
                                                        subarray.size(),
                                                        subarray[subarray.size()-1],
                                                        op,
                                                        identity_element);

           

            T processor_sum = 0;
           // #pragma omp parallel for reduction(+:processor_sum)
            for(int p=0; p<i; p++)
                processor_sum = op(processor_sum, sum[p]);
           
            
            if(i==0){
                    output[subarray_length*i] = subarray_scan[0];
                }
                else{
                    output[subarray_length*i] = op(subarray_scan[0], processor_sum);
                }
          // wrtiting final output corresponding to each subarray
            for(int k=1; k <subarray_length; k++){

                output[subarray_length*i + k] = subarray_scan[k];
                if(i==0){
                    output[subarray_length*i + k] = subarray_scan[k];
                }
                else{
                    output[subarray_length*i + k] = op(subarray_scan[k], processor_sum);
                }
            }

            }


   }
// converting prescan output to scan
    vector<T> scan = vector<T>(output.begin()+1,output.begin()+n);
    scan.push_back(op(scan[scan.size()-1],last_element));

// copying the scan to "out"
    copy(scan.begin(),scan.begin()+n_original,out);
    }

#endif // A0_HPP
 
