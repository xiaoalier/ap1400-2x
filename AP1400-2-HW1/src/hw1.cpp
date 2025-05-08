#include "hw1.h"
#include <vector>
#include <iostream>

namespace algebra{
    Matrix zeros(size_t n, size_t m){
        Matrix matrix;
        std::vector<double>sum(m,0.0);
        for (int i = 0; i < n; ++ i) matrix.push_back(sum);
        return matrix;
    }
    Matrix ones(size_t n, size_t m){
        Matrix matrix;
        std::vector<double>sum(m,1.0);
        for (int i = 0; i < n; ++ i) matrix.push_back(sum);
        return matrix;
    }
    Matrix random(size_t n, size_t m, double min, double max){
        if (min >= max) throw std::logic_error("min must small than max");
        Matrix matrix;
        std::random_device rd;  // 随机种子（硬件熵源）
        std::mt19937 gen(rd()); // 初始化随机引擎
        std::uniform_real_distribution<> dist(min, max); // 范围 [1, 100]
        int random_double = dist(gen);
        std::vector<double>sum(m,random_double);
        for (int i = 0; i < n; i ++) matrix.push_back(sum);
        return matrix;
    }
    void show(const Matrix& matrix){
        for (const auto& x : matrix){
            for (const auto& y : x){
                std::cout << std::fixed << std::setprecision(3)<<y<<" "; 
            }
            std::cout << '\n';
        }
    }
    Matrix multiply(const Matrix& matrix, double c){
        Matrix matrix_new;
        int m = matrix[0].size();
        for (const auto& x : matrix){
            std::vector<double>multiply(m,c);
            for (int i = 0; i < m; ++ i){
                multiply[i]*=x[i]; 
            }
            matrix_new.push_back(multiply);
        }
        return matrix_new;
    }
    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2){
        if (matrix1.size() == 0) return matrix2;
        else if (matrix2.size() == 0) return matrix1;
        int n1 = matrix1.size(), n2 = matrix2.size(), m1 = matrix1[0].size(), m2 = matrix2[0].size();
        if (m1 != n2) throw std::invalid_argument("this two matrix can not multiply");
        Matrix matrix = zeros(n1,m2);
        for (int i = 0; i < n1; ++ i){
            for (int j = 0; j < m2; ++ j){
                for (int k = 0; k < m1; k ++) matrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
        return matrix;
    }
    Matrix sum(const Matrix& matrix, double c){
        if (matrix.size() == 0) return matrix;
        Matrix matrix_new(matrix);
        for (auto& x : matrix_new){
            for (auto& y : x){
                y += c;
            }
        }
        return matrix_new;
    }
    Matrix sum(const Matrix& matrix1, const Matrix& matrix2){
        if (matrix1.size() == 0 and matrix2.size() == 0) return matrix1;
        else if (matrix1.size() == 0 or matrix2.size() == 0)throw std::logic_error("this two can not add");
        int n1 = matrix1.size(), n2 = matrix2.size(), m1 = matrix1[0].size(), m2 = matrix2[0].size();
        if (n1 != n2 or m1 != m2) throw std::logic_error("this two can not add");
        Matrix matrix(matrix1);
        for (int i = 0; i < n1; ++ i){
            for (int j = 0; j < m1; ++ j){
                matrix[i][j] += matrix2[i][j];
            }
        }
        return matrix;
    }
    Matrix transpose(const Matrix& matrix){
        if (matrix.size() == 0) return matrix;
        int n = matrix.size();
        int m = matrix[0].size();
        Matrix matrix_new = zeros(m,n);
        for (int i = 0; i < m ; ++ i){
            for (int j = 0; j < n; ++ j){
                matrix_new[i][j] = matrix[j][i];
            }
        }
        return matrix_new;
    }
    Matrix minor(const Matrix& matrix, size_t n, size_t m){
        int n1 = matrix.size();
        int m1 = matrix[0].size();
        Matrix matrix_new;
        for (int i = 0; i < n1; ++ i){
            if (i == n) continue;
            std::vector<double> minor;
            for (int j = 0; j < m1; ++ j){
                if (j == m) continue;
                minor.push_back(matrix[i][j]);
            }
            matrix_new.push_back(minor);
        }
        return matrix_new;
    }
    double determinant(const Matrix& matrix){
        if (matrix.size() == 0) return 1;
        if (matrix.size() != matrix[0].size()) throw std::logic_error("this matrix can not get determinat");
        int n = matrix.size();
        if (n == 1) return matrix[0][0];
        if (n == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        double ans = 0.0;
        for (int i = 0; i < n; ++ i){
            Matrix matrix_new = minor(matrix,0,i);
            ans += ((i % 2 == 0 ? 1 : -1) * matrix[0][i] * determinant(matrix_new));
        }
        return ans;
    }
    static Matrix adjoint(const Matrix& matrix){
        int n = matrix.size();
        Matrix matrix_new = zeros(n,n);
        for (int i = 0; i < n; ++ i){
            for (int j = 0; j < n; ++ j){
                int k = (i + j) & 1 ? -1 : 1;
                Matrix matrix_new2 = minor(matrix,i,j);
                matrix_new[i][j] = k * determinant(matrix_new2);
            }
        }
        return transpose(matrix_new);
    }
    Matrix inverse(const Matrix& matrix){
        if (matrix.size() == 0) return matrix;
        int n1 = matrix.size(), m1 = matrix[0].size();
        if (n1 != m1) throw std::logic_error("it must be suqare");
        if (determinant(matrix) == 0) throw std::logic_error("it's det must not be zero");
        Matrix matrix_new = adjoint(matrix);
        return multiply(matrix_new,1/determinant(matrix));
    }
    
    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis){
        int n1,m1,n2,m2;
        n1 = matrix1.size();
        n2 = matrix2.size();
        if (n1 == 0) return matrix2;
        if (n2 == 0) return matrix1;
        Matrix matrix_new;
        m1 = matrix1[0].size(), m2 = matrix2[0].size();
        if (axis == 0){
            if (m1 != m2) throw std::logic_error("wrong enter");
            matrix_new = zeros(n1 + n2, m1);
            for (int i = 0; i < n1 + n2; ++ i){
                for (int j = 0; j < m1; ++j){
                    matrix_new[i][j] = i < n1 ? matrix1[i][j] : matrix2[i - n1][j];
                }
            }
        }
        else if (axis == 1){
            if (n1 != n2) throw std::logic_error("wrong enter");
            matrix_new = zeros(n1, m1 + m2);
            for (int i = 0; i < n1; ++ i){
                for (int j = 0; j < m1 + m2; ++ j){
                    matrix_new[i][j] = j < m1 ? matrix1[i][j] : matrix2[i][j - m1];
                }
            }
        }
        else throw std::invalid_argument("axis must be 0 or 1");
        return matrix_new;
    }
    Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2){
        int n = matrix.size();
        if (n <= r1 or n <= r2) throw std::logic_error("wrong enter");
        if (n == 0) return matrix;
        Matrix matrix_new(matrix);
        std::swap(matrix_new[r1],matrix_new[r2]);
        return matrix_new;
    }
    Matrix ero_multiply(const Matrix& matrix, size_t r, double c){
        int n = matrix.size();
        if (n == 0) return matrix;
        if (n <= r) throw std::logic_error("wrong enter");
        Matrix matrix_new(matrix);
        for (auto& x : matrix_new[r]) x *=  c;
        return matrix_new;
    }
    Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2){
        int n = matrix.size();
        if (n <= r1 or n <= r2) throw std::logic_error("wrong enter");
        if (n == 0) return matrix;
        Matrix matrix_new(matrix);
        for (int i = 0; i < matrix[r2].size(); ++ i) matrix_new[r2][i] += matrix[r1][i] * c;
        return matrix_new;
    }
    Matrix upper_triangular(const Matrix& matrix){
        int n = matrix.size();
        if (n == 0) return matrix;
        
        int m = matrix[0].size();
        if (n != m) throw std::logic_error("mat must be square");
        if (n == 1) return matrix;
        Matrix matrix_new(matrix);
        for (int i = 0; i < n; ++ i){
            if (abs(matrix_new[i][i]) < 1e-9){
                bool swapped = false;
                for (int k = i + 1; k < n; ++ k){
                    if (abs(matrix_new[k][i]) > 1e-9){
                        swap(matrix_new[i],matrix_new[k]);
                        swapped =true;
                        break;
                    }
                }
                if (!swapped) continue;
            }
            double de = matrix_new[i][i];
            for (int j = i + 1; j < n; ++ j){
                double elem = matrix_new[j][i];
                matrix_new = ero_sum(matrix_new,i,-elem/de,j);
            }
        }
        return matrix_new;
    }
}