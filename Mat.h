#include<iostream>
#include <ctime>
#include<vector>
using namespace std;

class Matrx
{
public:
	Matrx(int rows, int cols);
	Matrx(int rows, int cols, float val);
	void divid(int, int, int, int, Matrx*);
	void put_val(int, int, int, int, Matrx&);
	~Matrx();
	int rows;
	int cols;
	vector<float>* mat_1 = new vector<float>;
	Matrx* Reg_Mul(Matrx&);
	void Mat_Equal(Matrx*);
	bool Check_Equality(Matrx&);
	Matrx* Mat_Add(Matrx*);
	Matrx* Strassen_Mul(Matrx*);


private:
	void Mat_Initialize();
	void Mat_Initialize(float);
};

Matrx::Matrx(int rows, int cols)
{
	this->rows = rows;
	this->cols = cols;

	Mat_Initialize();
}

Matrx::Matrx(int rows, int cols, float crt_val)
{
	this->rows = rows;
	this->cols = cols;

	Mat_Initialize(crt_val);
}

void Matrx::Mat_Initialize()
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			mat_1[i][j] = 0;
	}
}
void Matrx::Mat_Initialize(float crt_val)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			mat_1[i][j] = crt_val;
		}
	}
}


void Matrx::divid(int rs, int re, int cs, int ce, Matrx* mat2)
{
	for (int i = 0; rs < re; rs++, i++)
	{
		int _cs = cs;
		for (int j = 0; _cs < ce; _cs++, j++)
		{
			mat_1[i][j] = mat2->mat_1[rs][_cs];
		}
	}
}


void Matrx::put_val(int rs, int re, int cs, int ce, Matrx& mat2)
{
	for (int i = 0; rs < re; rs++, i++)
	{
		int _cs = cs;
		for (int j = 0; _cs < ce; _cs++, j++)
		{
			mat_1[rs][_cs] = mat2.mat_1[i][j];
		}
	}
}


void Matrx::Mat_Equal(Matrx* mat2)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			mat_1[i][j] = mat2->mat_1[i][j];
		}
	}
}

bool Matrx::Check_Equality(Matrx& mat2)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (mat_1[i][j] != mat2.mat_1[i][j])
				return false;
		}
	}
	return true;
}

Matrx* Matrx::Mat_Add(Matrx* mat2)
{
	if (rows != mat2->rows || cols != mat2->cols)
		cout << "Dimensions did not match" << endl;

	Matrx* result = new Matrx(rows, cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			result->mat_1[i][j] = mat_1[i][j] + mat2->mat_1[i][j];
		}
	}
	return result;
}



Matrx* Matrx::Reg_Mul(Matrx& mat2)
{
	if (cols != mat2.rows)
		cout << "Dimensions did not match" << endl;

	Matrx *result = new Matrx(rows, mat2.cols);

	for (int i = 0; i < result->rows; i++)
	{
		for (int j = 0; j < result->cols; j++)
		{
			result->mat_1[i][j] = 0;
			for (int k = 0; k < cols; k++)
			{
				result->mat_1[i][j] += mat_1[i][k] * mat2.mat_1[k][j];
			}
		}
	}

	return result;
}

Matrx* Matrx::Strassen_Mul(Matrx* mat2)
{


	if (cols != mat2->rows)
	{
		cout << "Dimensions did not match" << endl;
		return NULL;
	}

	if (rows == 2)
	{
		return this->Reg_Mul( *mat2);
	}
	else
	{
		int half = rows / 2;
		int end = rows;

		Matrx* result = new Matrx(rows, cols, 0);
		Matrx A11(half, half, 0);
		A11.divid(0, half, 0, half, this);
		Matrx A12(half, half, 0);
		A12.divid(0, half, half, end, this);
		Matrx A21(half, half, 0);
		A21.divid(half, end, 0, half, this);
		Matrx A22(half, half, 0);
		A22.divid(half, end, half, end, this);

		Matrx* B11 = new Matrx(half, half, 0);
		B11->divid(0, half, 0, half, mat2);
		Matrx* B12 = new Matrx(half, half, 0);
		B12->divid(0, half, half, end, mat2);
		Matrx* B21 = new Matrx(half, half, 0);
		B21->divid(half, end, 0, half, mat2);
		Matrx* B22 = new Matrx(half, half, 0);
		B22->divid(half, end, half, end, mat2);

		result->put_val(0, half, 0, half, *(A11.Strassen_Mul(B11)->Mat_Add( (A12.Strassen_Mul(B21)))));
		result->put_val(0, half, half, end, *(A11.Strassen_Mul(B12)->Mat_Add(A12.Strassen_Mul(B22))));
		result->put_val(half, end, 0, half, *(A21.Strassen_Mul(B11)->Mat_Add (A22.Strassen_Mul(B21))));
		result->put_val(half, end, half, end, *(A21.Strassen_Mul(B12)->Mat_Add (A22.Strassen_Mul(B22))));

		return result;
	}
}




Matrx::~Matrx()
{
	delete mat_1;
}


/*Testing Class*/
class Test
{
public:
	bool initialize_test();
	bool Reg_Mul_Test();
	bool Strassen_Mul_Test();
	bool Test::Equality_Test();
};
bool Test::initialize_test()
{
	bool result = true;
	Matrx first_mat(16, 16, 0);
	for (int i = 0; i < first_mat.rows; i++)
	{
		for (int j = 0; j < first_mat.cols; j++)
		{
			result &= (first_mat.mat_1[i][j] == 0);
		}
	}
	Matrx second_mat(16, 16, 2);
	for (int i = 0; i < second_mat.rows; i++)
	{
		for (int j = 0; j < second_mat.cols; j++)
		{
			result &= (second_mat.mat_1[i][j] == 2);
		}
	}
	return result;
}
bool Test::Equality_Test()
{
	bool temp = true;
	Matrx first_mat(16, 16, 1);
	Matrx second_mat(16, 16, 1);
	temp = first_mat.Check_Equality(second_mat);
	return temp;
}
bool Test::Reg_Mul_Test()
{
	bool temp = true;
	Matrx first_mat(16, 16);
	Matrx second_mat(16, 16);
	first_mat.Reg_Mul(second_mat);
	temp=first_mat.Check_Equality(second_mat);
	return temp;
}

bool Test::Strassen_Mul_Test()
{
	bool temp = true;
	Matrx first_mat(16, 16);
	Matrx second_mat(16, 16);
	Matrx* mat1 = first_mat.Reg_Mul(second_mat);
	Matrx* mat2 = first_mat.Strassen_Mul(&second_mat);
	temp &= (mat1->Check_Equality(*mat2));
return temp;
}
