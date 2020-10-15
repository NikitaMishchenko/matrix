#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <iostream>

class matrix
{
private:
    static char out_separator;
protected:
    int height;
    int width;
    double** data; ///����� ��������� ������ �� vector
public:

///constructors//des//copy
    matrix(); ///empty
    void matrix_alloc(int , int); ///�������� ������ �������� �������
    matrix(int , int); ///rectangle
    matrix(int); ///Square matrix
    matrix(int n, char);///���������� 'E' - ���������, 'n' - � �������� ����������� ����� � 1 �� n.
    matrix(int i, int j, char);///'E' - ���������, 'n' - � �������� ����������� ����� � 1 �� n.
    ~matrix(); ///destructor
    matrix(const matrix&); ///copy

///primary
    void delete_data();///������� ������ ������� ������� ������ � �������� 0�0 /// ����� ��������������
    int get_height() const;
    int get_width() const;
    double get_data(int i, int j) const;
    void set_height(int h);
    void set_width(int w);
    void set_data(int i, int j, double d);

    void set_out_separator(char s);

    void set_column_data(int h, int w, double data);///��������� ��������� ������� � 0 �� h, ������� � �������� w, ���� ��� ������ - �������.

    void make_matrix_data(int h, int w, double d);

    void insert_array_to_matrix_column(int index_column, double* arr, int length);///overwrite and fit to matrix.height();
    void insert_array_to_matrix_row(int index_row, double* arr, int length);///overwrite and fit to matrix.width();

    void get_array_from_column(int index_column, double* arr);///length_of_array = matrix.height
    void get_array_from_row(int index_row, double* arr);///length_of_array = matrix.width

    double get_max_el_column(int index_column);
        double get_max_el_column_gap(int index_column, int gap_index_from, int gap_index_to);
    double get_min_el_column(int index_row);

    int get_max_elind_column_gap(int index_column, int gap_index_from, int gap_index_to);
    int get_min_elind_column_gap(int index_column, int gap_index_from, int gap_index_to);

///operatros
    matrix& operator=(const matrix& A);

    matrix& operator+=(const matrix& A);
    friend matrix operator+(const matrix& A, const matrix& B);

    matrix operator* (double a); /// ����� �������������� /// ����������� �/� �������� ��������� O(n^2.81)
    friend matrix operator*(const double& a, const matrix& A);
    friend matrix operator*(const matrix& A, const matrix& B);
    matrix operator*=(const matrix& A);

    matrix operator/ (const double a);

    matrix operator-= (const matrix& A);
    friend matrix operator-(const matrix& A, const matrix& B);

    friend std::ostream& operator << (std::ostream &C, const matrix &A);
    friend std::istream& operator >> (std::istream &C, matrix &A);

///FUNCTIONS
///math
    matrix abs_column(int index_column);
    void abs_index_column(int index_column);
    //matrix abs_row( int index_row);
    //matrix abs();///all

    matrix transpose();
        void transpose(matrix& ); /// ��� ��������� � ����� �������?
    matrix minor(int l,int k); ///"����������" ����� ������� � ������ /// ���������� �������!!!!!
    double determinant(); ///default algorithm  /// minor /// ����� ��������������
        void determinant(double& determinant);
    matrix inverse();
        void inverse_this();
    matrix elementwise_division(matrix B);/// sizes1==sizes2 /// A/B == A.elementwise_division(B)

    int width_coord_of_max_elem_from_row(int column_index);///coordinates counts from 0 //?
    void triangulation(); ///?
    void triangulation1(matrix& B); ///elem err ///?

    matrix average_column();

    void column_multiply(int index_column, double multiplier);
    void column_shift(int index_column, double sh);///��������� � ������� �������� ������� index_column sh
        matrix m_column_shift(int index_column, double sh);///����������� �������

    matrix matrix_exponential(int index_order);
    void Fibonachi2x2(int f_number);

    matrix sort_by_column_index(int index_column);

///structure
    void erase_row(int h);/// - /// ������� ������ � ������ ������� /// ����� ��������������
    void erase_column(int w);///| /// ������� ������� � ������ ������� /// ����� ��������������
    matrix get_row(int h);/// ���������� ������ ������ ������� �� ������ /// ����� ��������������
    matrix get_column(int w);
    matrix get_matrix_part(int h1, int h2, int w1, int w2); ///simple part from h1,w1 to h2,w2 ///���������� ����� ������� ������� �� ������� h1 �� h2 � �� ������ w1 �� w2

    matrix merge_height(const matrix B, int pozition_row); ///��������
        matrix merge_height(matrix B);///� �����
    matrix merge_width(matrix B, int pozition_column); ///��������
        matrix merge_width(matrix B);///� �����
        void merge_height(matrix A, matrix B, int pozition_row); /// �� ��������!!!
        void merge_width(matrix A, matrix B, int pozition_column); //

    void swap_elements(int i1, int j1, int i2, int j2);///����� ������� ����� ��������� /// ����� ��������������
        void swap_columns(int w1, int w2);
        void swap_rows(int h1, int h2);

///types
    //matrix array_to_matrix_column(double* arr, int length);
   // matrix array_to_matrix_row(double* arr, int length);
    void array_to_matrix_column(double* arr, int length);
    matrix matrix_column_from_array(double* arr, int length);
    void array_to_matrix_row(double* arr, int length);
    matrix matrix_row_from_array(double* arr, int length);

///files
    void Find_Size_of_matrix_file(const std::string file_name); ///������������� height � width ������ ������� �� �������� ����� � �������, ����� ����� �� ����� /// ����� ��������������
    void Find_Size_of_matrix_emptstringend_file(const std::string file_name, int emptstringend_num);///����� � ����� ����� ������ ������
    void Load_matrix(const std::string file_name); ///memory allocation
    void Load_matrix_emptstringend(const std::string file_name, int emptstringend_num);
        //void Load_matrix_with_header(const string file_name);///������� /// ������� � ��������� ����
    void Save_matrix(const std::string file_name);
    void Save_matrix_end(const std::string file_name); /// ���������� ������� � ����� �����

///extra
    void info();
};

///outer
    bool FileExists(std::string fname);
    bool FileExists1(std::string fname);
    bool FileExists2(const char *fname);

#endif // MATRIX_H_INCLUDED
