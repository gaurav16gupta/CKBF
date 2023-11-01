#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bitArray.h"

using namespace std;

bitArray::bitArray(int size){
  ar_size = size;
  A = new char[ar_size/8 +1];
  for (int i=0; i<(ar_size/8 +1); i++ ){
    A[i] = '\0'; // Clear the bit array
    }
}

void bitArray::clearall(void) {
  for (int i=0; i<(ar_size/8 +1); i++ ){
    A[i] = '\0'; // Clear the bit array
    }
}
  
void bitArray::SetBit(uint k) {
   A[(k/8)] |= (1 << (k%8));
}

void bitArray::ClearBit(uint k) {
   A[(k/8)] &= ~(1 << (k%8));
}

bool bitArray::TestBit(uint k) {
   return (A[(k/8)] & (1 << (k%8)));
}

void bitArray::ANDop(char* B){
  for (int i=0; i<(ar_size/8 +1); i++ ){
    A[i] &= B[i];
  }
}

int bitArray::getcount(void){
  int count = 0;

  for (int kp=0; kp<ar_size; kp++ ){
    // std::cout << "kp is" <<kp<< ' ';
    if (A[(kp/8)] & (1 << (kp%8)) ){
      count++;
    }
  }
  return count;
}

vector<uint> bitArray::get1locs(void){
  vector<uint> locs;

  for (int kp=0; kp<ar_size; kp++ ){
    // std::cout << "kp is" <<kp<< ' ';
    if (A[(kp/8)] & (1 << (kp%8)) ){
      locs.push_back(kp);
    }
  }
  return locs;
}

void bitArray::serializeBitAr(string BF_file){
  ofstream out;
  out.open(BF_file);

  if(! out){
    cout<<"Cannot open output file\n";
  }
  out.write(A,ar_size/8 +1);
    out.close();
}

void bitArray::deserializeBitAr(std::vector<string> BF_file){
  ifstream in(BF_file[0]);
  // cout<<BF_file[0]<<endl;
  if(! in){
    cout<<"Cannot open input file\n";
  }
  in.read(A,ar_size/8 +1); //optimise it
  in.close();
}
   // // Check if SetBit() works:
   //
   // for ( i = 0; i < 320; i++ )
   //    if ( TestBit(A, i) )
   //       printf("Bit %d was set !\n", i);
   //
   // printf("\nClear bit poistions 200 \n");
   // ClearBit( A, 200 );
   //
   // // Check if ClearBit() works:
   //
   // for ( i = 0; i < 320; i++ )
   //    if ( TestBit(A, i) )
   //       printf("Bit %d was set !\n", i);
