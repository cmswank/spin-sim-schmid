#include "BField.h"
#include "BList.h"

#include <iostream>

using namespace std;

int test_BFieldAdder();
int test_BFieldAdder2();
void test_BFieldScaler();

int main(int argc, char **argv) {
  test_BFieldAdder2();
  test_BFieldScaler();
}

void test_BFieldScaler() {
  double temp_field1[] = {10., -1, 2};
  double f = .125;
  BFieldConst *field1 = new BFieldConst(temp_field1);
  BFactor *factor = new BFactor(f);
  BField *field = new BFieldScaler(field1, factor);
  double temp[3];
  BFieldVars vars;
  field->getField(temp, vars);
  f = factor->getFactor(vars);
  delete field;
  for (int i=0; i<3; i++) {
    temp_field1[i] *= f;
    double diff = temp[i] - temp_field1[i];
    if (diff != 0) {
      cout << "!!! BFieldScaler Failed" << endl;
      cout << "\t" << temp[i] << " - " << temp_field1[i]
	   << " = " << diff << endl;
      return;
    }
  }
  cout << "OK BFieldScaler" << endl;
}
int test_BFieldAdder2() {
  double temp_field1[] = {10, -1, 2};
  double temp_field2[] = {-5, 6, 3};
  double result[] = {10, 10, 10};
  double temp[3];
  BFieldConst *field1 = new BFieldConst(temp_field1);
  BFieldConst *field2 = new BFieldConst(temp_field2);
  BFieldConst *field1B = new BFieldConst(temp_field1);
  BFieldConst *field2B = new BFieldConst(temp_field2);
  BFieldAdder *adder = new BFieldAdder();
  adder->append(field1);
  adder->append(field2);
  BFieldAdder *adder2 = new BFieldAdder();
  adder2->append(field1B);
  adder2->append(field2B);
  BField *field = new BFieldAdder(adder, adder2);
  BFieldVars vars;
  field->getField(temp, vars);
  delete field;
  for (int i=0; i<3; i++) {
    if (temp[i] != result[i]) {
      cout << "!!! BFieldAdder test did not succeed" << endl;
      return 0;
    }
  }
  cout << "OK BFieldAdder" << endl;
  return 1;
}

int test_BFieldAdder() {
  double temp_field1[] = {10, -1, 2};
  double temp_field2[] = {-5, 6, 3};
  double result[] = {10, 10, 10};
  double temp[3];
  BField *field1 = new BFieldConst(temp_field1);
  BField *field2 = new BFieldConst(temp_field2);
  BFieldAdder *adder = new BFieldAdder();
  adder->append(field1);
  adder->append(field2);
  BFieldAdder *adder2 = new BFieldAdder();
  adder2->append(field1);
  adder2->append(field2);
  BField *field = new BFieldAdder(adder, adder2);
  BFieldVars vars;
  field->getField(temp, vars);
  for (int i=0; i<3; i++) {
    if (temp[i] != result[i]) {
      cout << "!!! BFieldAdder test did not succeed" << endl;
      return 0;
    }
  }
  cout << "OK BFieldAdder" << endl;
  return 1;
}
