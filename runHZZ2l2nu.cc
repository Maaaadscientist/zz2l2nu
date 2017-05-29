//-*- c-basic-offset: 4; -*-
#include <iostream>
#include <thread>
#include <TString.h>
#include <TChain.h>
#include "Includes/HZZ2l2nuLooper.h"

using namespace std;

int main(int argc, char **argv)
{
  HZZ2l2nuLooper myHZZlooper("mySkim.root"); //please add the parameter that you think are needed
  myHZZlooper.Loop(); //same here 
}
