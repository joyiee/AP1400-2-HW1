#include "hw1.h"

void valid_check(size_t n, size_t m)
{
    if(!n || !m)
    {
        throw std::logic_error("Invalid matrix size!");
    }
}

void valid_check(const Matrix &matrix)
{
    size_t column = matrix[0].size();
    for(auto row: matrix)
    {
        if(row.size()!= column)
           throw std::logic_error("Invalid matrix");
    }
    return;
}

bool valid_check(const Matrix &matrix1, const Matrix &matrix2)
{
    bool c1 = matrix1.empty(), c2 = matrix2.empty();
    if((c1&&!c2)||(!c1&&c2))
        throw std::logic_error("Invalid matrix");
    return (c1&&c2);
}

double generateRandomNumber(double min, double max)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(min, max);

    return dis(gen);
}

namespace algebra{
    Matrix zeros(size_t n, size_t m){
        valid_check(n, m);
        Matrix res(n, std::vector<double> (m, 0));
        return res;
    }

    Matrix ones(size_t n, size_t m){
        valid_check(n, m);
        Matrix res(n, std::vector<double> (m, 1));
        return res;
    }

    Matrix random(size_t n, size_t m, double min, double max)
    {
        valid_check(n, m);
        if(min > max){
            throw std::logic_error("Invalid bound");
        }

        Matrix res(n, std::vector<double> (m, 0));
        for(size_t i = 0; i < n; i++)
        {
            for(size_t j = 0; j < m; j++)
            {
                res[i][j] = generateRandomNumber(min, max);
            }
        }
        return res;
    }

    void show(const Matrix& matrix)
    {
        if(matrix.empty())
            return;
        valid_check(matrix);
        std::cout << std::showpoint << std::setprecision(3);
        for(auto &row: matrix)
        {
            for(auto &elem: row)
            {
                std::cout << std::setw(9) << elem << " ";
            }
            std::cout << "\n";
        }
    }

    Matrix multiply(const Matrix& matrix, double c)
    {
        if(matrix.empty())
        {
            Matrix res;
            return res;
        }
        valid_check(matrix);
        Matrix res = matrix;
        for(auto &row: res)
        {
            for(auto &elem : row)
            {
                elem *= c;
            }
        }
        return res;
    }

    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2)
    {
        if(valid_check(matrix1, matrix2)){
            Matrix res;
            return res;
        }
        valid_check(matrix1);
        valid_check(matrix2);
        
        if(matrix1[0].size() != matrix2.size())
            throw std::logic_error("Invalid matrix size when multiplying");
        
        Matrix res;
        
        for(size_t i = 0; i < matrix1.size(); i++)
        {
            std::vector<double> line;
            for(size_t j =0; j < matrix2[0].size(); j++)
            {
                double sum = 0;
                for(size_t k = 0; k < matrix1[0].size(); k++)
                {
                    sum+= matrix1[i][k] * matrix2[k][j];
                }
                line.push_back(sum);
            }
            res.push_back(line);
        }
        return res;
    }

    Matrix sum(const Matrix& matrix, double c)
    {
        if(matrix.empty())
        {
            Matrix res;
            return res;
        }
        valid_check(matrix);
        Matrix res = matrix;
        for(auto &row : res)
        {
            for(auto &elem : row)
            {
                elem += c;
            }
        }
        return res;
    }

    Matrix sum(const Matrix& matrix1, const Matrix& matrix2)
    {
        if (valid_check(matrix1, matrix2))
        {
            Matrix res;
            return res;
        }
        valid_check(matrix1);
        valid_check(matrix2);
        int r1 = matrix1.size(), r2 = matrix2.size();
        int c1 = matrix1[0].size(), c2 = matrix2[0].size();
        if(!(r1==r2 && c1 == c2))
        {
            throw std::logic_error("Invalid matrix size when adding up two matrixs");
        }
        Matrix res = matrix2;
        for(int i = 0; i < r1; i++)
        {
            for(int j = 0; j < c1; j++)
            {
                res[i][j]+= matrix1[i][j];
            }
        }
        return res;
    }
    
    Matrix transpose(const Matrix & matrix)
    {
        if(matrix.empty())
        {
            Matrix res;
            return res;
        }
        valid_check(matrix);

        int row = matrix.size(), column = matrix[0].size();
        Matrix res = zeros(column, row);
        for(int i = 0; i < column ; i++)
        {
            for(int j = 0; j < row; j++)
            {
                res[i][j] = matrix[j][i];
            }
        }
        return res;
    }

    Matrix minor(const Matrix& matrix, size_t n, size_t m)
    {
        valid_check(matrix);
        int row = matrix.size(), column = matrix[0].size();
        if(!(n < row && m < column))
            throw std::logic_error("Invalid error");
        Matrix res;
        for(int i = 0; i < row ; i++)
        {
            if(i == n) continue;
            std::vector<double> line;
            for(int j = 0; j< column; j++)
            {
                if(j == m) continue;
                line.push_back(matrix[i][j]);
            }
            res.push_back(line);
        }
        return res;
    }

    double determinant(const Matrix& matrix)
    {
        if(matrix.empty())
            return 1;
        valid_check(matrix);
        int row = matrix.size(), column = matrix[0].size();
        if(row != column)
            throw std::logic_error("Invalid matrix when calculating determinant");
        if(row == 1)
            return matrix[0][0];
        double op = 1, res = 0;
        for(int i = 0; i < column; i++)
        {
            res += op*matrix[0][i]*determinant(minor(matrix, 0, i));
            op *= -1;
        }
        return res;
    }

    Matrix inverse(const Matrix& matrix)
    {
        if(matrix.empty())
        {
            Matrix res;
            return res;
        }
        double D = determinant(matrix);
        if(D == 0)
            throw std::logic_error("Unable to get inverse matrix");
        
        int n = matrix.size();
        Matrix res = zeros(n, n);
        
        for(size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                res[i][j] = determinant(minor(matrix, j, i)) / D;
                if((i+j)%2)
                {
                    res[i][j] *= -1;
                }
            }
        }  
        return res;
    }
    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis)
    {
        if(axis != 0 && axis!= 1)
        {
            throw std::logic_error("Invalid axis");
        }
        if(valid_check(matrix1, matrix2))
        {
            Matrix res;
            return res;
        }
        valid_check(matrix1);
        valid_check(matrix2);

        int r1 = matrix1.size(), r2 = matrix2.size();
        int c1 = matrix1[0].size(), c2 = matrix2[0].size();
        
        if((axis == 0 && c1!= c2)||(axis==1 && r1!=r2))
        {
            throw std::logic_error("Invalid size when concatenating");
        }

        Matrix res;

        if(axis == 0)
        {
            for(auto line : matrix1)
            {
                res.push_back(line);
            }
            for(auto line : matrix2)
            {
                res.push_back(line);
            }
        }
        else
        {
            for(size_t i = 0; i< r1; i++)
            {
                std::vector<double> line = matrix1[i];
                for(auto elem : matrix2[i])
                {
                    line.push_back(elem);
                }
                res.push_back(line);
            }
        }
        return res;
    }

    Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2)
    {
        if(matrix.empty())
        {
            Matrix res;
            return res;
        }
        valid_check(matrix);
        size_t row = matrix.size();
        if(r1 < 0 || r2 < 0 || r1 >= row || r2 >=row)
            throw std::logic_error("Out of bound");
        Matrix res = matrix;
        swap(res[r1], res[r2]);
        return res;
    }

    Matrix ero_multiply(const Matrix& matrix, size_t r, double c)
    {
        if(matrix.empty())
        {
            Matrix res;
            return res;
        }
        valid_check(matrix);
        size_t row = matrix.size();
        if(r < 0 || r>=row)
            throw std::logic_error("Out of bound");
        Matrix res = matrix;
        for(auto &elem : res[r])
        {
            elem *= c;
        }
        return res;
    }
    Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2)
    {
        if(matrix.empty())
        {
            Matrix res;
            return res;
        }
        valid_check(matrix);
        size_t row = matrix.size();
        if(r1 < 0 || r2 < 0 || r1 >= row || r2 >=row)
            throw std::logic_error("Out of bound");
        Matrix res = matrix;
        for (size_t i = 0; i < row; ++ i) {
            res[r2][i] += c * res[r1][i];
        }
        return res;
    }
    Matrix upper_triangular(const Matrix& matrix)
    {
        if(matrix.empty())
        {
            Matrix res;
            return res;
        }
        valid_check(matrix);
        int row = matrix.size(), column = matrix[0].size();
        if (row != column) {
            throw std::logic_error("Invalid size when upper_triangularing");
        }
        Matrix res = matrix;

        for (size_t i = 0; i < row - 1; ++ i) {
            for (size_t nonZero = i; nonZero < row; ++ nonZero) {
                if (res[nonZero][i]) {
                    swap(res[i], res[nonZero]);
                    break;
                }
            }
            for (size_t j = i + 1; j < row; ++ j) {
                double r = res[j][i] / res[i][i];
                for (size_t k = i; k < row; ++ k) {
                    res[j][k] -= r * res[i][k];
                }
            }
        }
        return res;
    }
}