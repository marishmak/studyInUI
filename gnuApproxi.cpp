//Mariia Shmakova
//m.shmakova@innopolis.university

#include <iostream>
#include <utility>
#include <vector>
#include <iomanip>
#include <cmath>
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"

using namespace std;

template <typename T>
class ColumnVector{
public:
    int n{};
    int m{};
    vector<vector<T>> matrix;

    ColumnVector() = default;

    ColumnVector(int n, int m, vector<vector<T>> matrix1) : n{n}, m{m}, matrix{std::move(matrix1)}
    {}

    ColumnVector(int n, vector<double> b) :  n{n}, m{1}{

        for(int i = 0; i < this->n; ++i) {
            vector<double> vector1(this->m, 0);
            this->matrix.push_back(vector1);
        }

        for(int i = 0; i < this->n; ++i) {
            for(int j = 0; j < this->m; ++j) {
                this->matrix[i][j] = b[i];
            }
        }

    }

    friend istream& operator >> (istream& input, ColumnVector &matrix1){

        input >> matrix1.n;
        matrix1.m = 1;

        vector<vector<T>> vector1(matrix1.n, vector<T>(matrix1.m, 0.00));
        matrix1.matrix = vector1;

        for(int i = 0; i < matrix1.n; ++i) {
            for(int j = 0; j < matrix1.m; ++j) {
                input >> matrix1.matrix[i][j];
            }
        }
        return input;
    }


    friend ostream& operator << (ostream& output, const ColumnVector &matrix1){

        for( int i=0; i<matrix1.n; ++i){
            output << fixed << setprecision(4);
            if(0.0001>abs(matrix1.matrix[i][0])){
                output << 0.00;
            }else{
                output << matrix1.matrix[i][0];
            }
            output << "\n";
        }
        return output;
    }

    ColumnVector& operator= (const ColumnVector& matrix1)= default;

};


template <typename T>
class Matrix {
public:
    int m{};
    int n{};
    vector<vector<T>>matrix;

    Matrix() = default;

    Matrix(int n, int m):n{n}, m{m}{
        vector<vector<T>> vector1(n, vector(m, 0.0));
        this->matrix = vector1;
    }

    Matrix(int n, int m, vector<vector<T>> matrix1) : n{n}, m{m}, matrix{std::move(matrix1)}
    {}

    Matrix(const Matrix &matrix1, const Matrix &matrix2) {
        this->n = matrix1.n;
        this->m = matrix1.m + matrix2.m;

        vector<vector<T>> augmented;
        for(int i = 0; i < n; ++i) {
            vector<T> aug(m, 0.00);
            augmented.emplace_back(aug);
        }

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < matrix1.m; ++j) {
                augmented[i][j] = matrix1.matrix[i][j];
            }
            for(int j = 0; j < matrix2.m; ++j) {
                augmented[i][j + matrix1.m] = matrix2.matrix[i][j];
            }
        }

        this->matrix = augmented;
    }

    friend istream& operator >> (istream& input, Matrix &matrix1){

        input >> matrix1.n;
        matrix1.m = matrix1.n;

        vector<vector<T>> vector1(matrix1.n, vector<T>(matrix1.m, 0.00));
        matrix1.matrix = vector1;


        for(int i = 0; i < matrix1.n; ++i) {
            for(int j = 0; j < matrix1.m; ++j) {
                input >> matrix1.matrix[i][j];
            }
        }
        return input;
    }


    friend ostream& operator << (ostream& output, const Matrix &matrix1){

        for( int i=0; i<matrix1.n; ++i){
            output << fixed << setprecision(4);

            if(0.0001>abs(matrix1.matrix[i][0])){
                output << 0.00;
            }else{
                output << matrix1.matrix[i][0];
            }

            for(int j=1; j<matrix1.m; ++j){
                output << fixed << setprecision(4);
                if(0.0001>abs(matrix1.matrix[i][j])){
                    output << " " << 0.00;
                }else{
                    output << " " << matrix1.matrix[i][j];
                }
            }
            output << "\n";
        }
        return output;
    }

    Matrix& operator= (const Matrix& matrix1)= default;

    auto operator* (const ColumnVector<T> &matrix1) {
        vector<vector<double> > mulMatrix;
        for(int i = 0; i < n; ++i) {
            vector<double> vector1(matrix1.m, 0.00);
            mulMatrix.emplace_back(vector1);
        }
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < matrix1.m; ++j) {
                for(int k = 0; k < m; ++k) {
                    mulMatrix[i][j] += (this->matrix[i][k] * matrix1.matrix[k][j]);
                }
            }
        }
        return *new ColumnVector(n, matrix1.m, mulMatrix);
    }

    auto operator* (const Matrix &matrix1) {
        vector<vector<T> > multiply(matrix1.n, vector<T>(matrix1.m, 0.00));

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < matrix1.m; ++j) {
                for(int k = 0; k < m; ++k) {
                    multiply[i][j] += (this->matrix[i][k] * matrix1.matrix[k][j]);
                }
            }
        }
        return *new Matrix(n, matrix1.m, multiply);
    }

    Matrix normalization () {
        vector<vector<T>> diagonal;

        for(int i = 0; i < m; ++i) {
            vector<T> diag(n, 0.00);
            diagonal.template emplace_back(diag);
        }
        for(int i = 0; i < n; ++i) {
            diagonal[i][i] = 1 / this->matrix[i][i];
        }

        return *new Matrix(n, m, diagonal);
    }

    auto transpose () {
        vector<vector<T> > transposeMatrix(m, vector<T>(n, 0.0));

        for(int i = 0; i < m; ++i) {
            for(int j = 0; j < n; ++j) {
                transposeMatrix[i][j] = this->matrix[j][i];
            }
        }
        return *new Matrix(m, n, transposeMatrix);
    }

};

template <typename  T>
class IdentityMatrix : public Matrix<T>{
public:
    int n;

    explicit IdentityMatrix(int n) : n{n}{
        vector<vector<T>> identity(n, vector<T> (n, 0));

        for (int i=0; i<n; ++i){
            for(int j=0; j<n; ++j){
                if(i==j){
                    identity[i][j] = 1;
                }else{
                    identity[i][j] = 0;
                }
            }
        }

        Matrix<T>::matrix = identity;
        Matrix<T>::m = n;
        Matrix<T>::n = n;
    }

    ~IdentityMatrix() = default;

    IdentityMatrix& operator= (const Matrix<T>& matrix1){
        Matrix<T>::matrix = matrix1.matrix;
        Matrix<T>::m = matrix1.n;
        Matrix<T>::n = matrix1.n;
        return *this;
    }
};

template <typename T>
class EliminationMatrix: public Matrix<T> {
public:
    int n{};

    EliminationMatrix() = default;
    EliminationMatrix(Matrix<T>& matrix, int n, int m) : n{matrix.n}{
        n++; m++;
        vector<vector<T>> elimination(matrix.n, vector<T> (matrix.n, 0.00));

        for (int i=0; i<this->n; ++i){
            for(int j=0; j<this->n; ++j){
                if(i==j){
                    elimination[i][j] = 1;
                }else{
                    elimination[i][j] = 0;
                }
            }
        }

        elimination[n-1][m-1] = -matrix.matrix[n-1][m-1];
        elimination[n-1][m-1] /= matrix.matrix[m-1][m-1];

        Matrix<T>::matrix = elimination;
        Matrix<T>::m = matrix.m;
        Matrix<T>::n = matrix.n;
    }

};

template<typename T>
class PermutationMatrix: public Matrix<T> {
public:

    int n{};


    PermutationMatrix() = default;

    PermutationMatrix(Matrix<T>& matrix, int n, int m) : n{matrix.n}{
        n++; m++;
        vector<vector<T>> permutation(matrix.n, vector<T> (matrix.n, 0.00));

        for (int i=0; i<this->n; ++i){
            for(int j=0; j<this->n; ++j){
                if(i==j){
                    permutation[i][j] = 1;
                }else{
                    permutation[i][j] = 0;
                }
            }
        }
        swap(permutation[n-1], permutation[m-1]);

        Matrix<T>::matrix = permutation;
        Matrix<T>::m = matrix.m;
        Matrix<T>::n = matrix.n;
    }

    ~PermutationMatrix() = default;

};


template <typename T>
void solvingEquation(Matrix<T>& A, ColumnVector<T>& CV) {
    int stepCounter = -1;
    T maxElement;
    bool mark;
    int maxElementIndex;
    ++stepCounter;

    for (int i=0; i<A.n; ++i){
        maxElement = abs(A.matrix[i][i]);
        maxElementIndex = -10^10;
        mark = false;
        for(int j=i; j<A.n; ++j){
            if(maxElement < abs(A.matrix[j][i])){
                maxElementIndex = j;
                maxElement =  abs(A.matrix[j][i]);
                mark = true;
            }
        }

        if(mark){
            auto Permutation = *new PermutationMatrix(A, i, maxElementIndex);
            auto permutatedA = Permutation * A;
            auto permutatedCV = Permutation * CV;
            A = permutatedA;
            CV = permutatedCV;
            ++stepCounter;
        }

        for(int j = i + 1; j < A.n; ++j) {
            if(0.0001<abs(A.matrix[j][i])) {
                auto Elimination = *new EliminationMatrix(A, j, i);
                auto eliminatedA = Elimination * A;
                auto eliminatedCV = Elimination * CV;
                A = eliminatedA;
                CV = eliminatedCV;
                ++stepCounter;
            }
        }
    }

    for(int i = A.n - 1; i > -1; --i) {
        for(int j = i - 1; j > -1; --j) {
            if(0.0001<abs(A.matrix[j][i])) {
                auto Elimination = *new EliminationMatrix(A, j, i);
                auto eliminatedA = Elimination * A;
                auto eliminatedCV = Elimination * CV;
                A = eliminatedA;
                CV = eliminatedCV;
                ++stepCounter;
            }
        }
    }


    auto Norm = A.normalization();
    auto normalizedA = Norm * A;
    auto normalizedCV = Norm * CV;
    A = normalizedA;
    CV = normalizedCV;
}

template <typename T>
auto inverse(Matrix<T> A){
    auto identity = *new IdentityMatrix<T>(A.n);
    int stepCounter = -1;
    T maxi;
    int indexMaxi = -10^10;
    bool mark = false;

    auto augmented = *new Matrix<T>(A, identity);
    ++stepCounter;


    for (int i=0; i<A.n; ++i){

        mark = false;
        maxi = A.matrix[i][i];
        for(int j=i+1; j<A.n; ++j){
            if(maxi < abs(A.matrix[j][i])){
                indexMaxi = j;
                maxi =  abs(A.matrix[j][i]);
                mark = true;
            }
        }

        if(mark){
            auto P = *new PermutationMatrix(A, i, indexMaxi);
            auto permutatedA = P * A;
            auto permutatedIdentity = P * identity;
            A = permutatedA;
            identity = permutatedIdentity;
            auto augmentedAndPermutated = *new Matrix(permutatedA, permutatedIdentity);
            ++stepCounter;
        }

        for(int j = i + 1; j < A.n; ++j) {
            if(A.matrix[j][i] != 0) {
                auto E = *new EliminationMatrix(A, j, i);
                auto eliminatedA = E * A;
                auto eliminatedIdentity = E * identity;
                A = eliminatedA;
                identity = eliminatedIdentity;
                auto augmentedAndEliminated = *new Matrix(eliminatedA, eliminatedIdentity);
                ++stepCounter;
            }
        }
    }

    for(int i = A.n - 1; i > -1; --i) {

        for(int j = i - 1; j > -1; --j) {
            if(A.matrix[j][i] != 0) {
                auto E = *new EliminationMatrix(A, j, i);
                auto eliminatedA = E * A;
                auto eliminatedIdentity = E * identity;
                A = eliminatedA;
                identity = eliminatedIdentity;
                auto augmentedAndEliminated = *new Matrix(eliminatedA, eliminatedIdentity);
                stepCounter++;
            }
        }
    }


    auto normalized = A.normalization();
    auto normalizedA = normalized * A;
    auto normalizedIdentity = normalized * identity;
    identity = normalizedIdentity;
    auto augmentedAndDiagonal = *new Matrix<T>(normalizedA, normalizedIdentity);

    return identity;
}


int main(){

    FILE* pipe = _popen(GNUPLOT_NAME, "w");
    int m;
    vector<double> t;
    vector<double> b;
    int n;
    double input;

    cin >> m;

    for(int i=0; i<m; ++i){
        cin >> input;
        t.emplace_back(input);
        cin >> input;
        b.emplace_back(input);
    }

    auto bColumned = *new ColumnVector<double>(b.size(), b);

    cin >> n;

    auto A = *new Matrix<double>(m, n+1);

    for(int i = 0; i < A.m; ++i) {
        for(int j = 0; j < A.n; ++j) {
            A.matrix[j][i] = pow(t[j], i);
        }
    }


    Matrix<double> transpose = A.transpose();
    Matrix<double> transposeA = transpose*A;
    IdentityMatrix<double> transposeAInverse = inverse(transposeA);
    ColumnVector<double> transposeB = transpose*bColumned;
    ColumnVector<double> answer = transposeAInverse*transposeB;


    cout << "A:\n" << A;
    cout << "A_T*A:\n" << transposeA;
    cout << "(A_T*A)^-1:\n" << transposeAInverse;
    cout << "A_T*b:\n" << transposeB;
    cout << "x~:\n" << answer;


    if(pipe != nullptr) {

        fprintf(pipe, "set terminal wxt size 800,600\n");
        fprintf(pipe, "plot [0 : 40] [0 : 200] ");

        for(int i = 0; i < answer.n; ++i) {
            if(i==0){
                fprintf(pipe, "%lf*x**%d", answer.matrix[0][0], 0);
            }else {
                fprintf(pipe, " + %lf*x**%d", answer.matrix[i][0], i);
            }
        }


        fprintf(pipe, " , '-' using 1:2 title 'given points' with points\n");
        for(int i = 0; i < m; ++i) {
            fprintf(pipe, "%f\t%f\n", t[i], b[i]);
            cout << t[i] << " " << b[i] << "\n";
        }

        fprintf(pipe, "e\n");
        fflush(pipe);
        _pclose(pipe);
    }

    return 0;

}

