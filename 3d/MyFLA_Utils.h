#include "FLAME.h"

int MyFLA_Obj_set_to_zero( FLA_Obj A );

int MyFLA_Obj_set_to_one( FLA_Obj A );

int MyFLA_Triu( FLA_Obj A, int num_diags );

int MyFLA_Tril( FLA_Obj A, int num_diags );

int MyFLA_Set_to_identity( FLA_Obj A );

int MyFLA_Set_main_diagonal_to_one( FLA_Obj A );

int MyFLA_Zero_strict_lower_triangular( FLA_Obj A );

int MyFLA_Zero_strict_upper_triangular( FLA_Obj A );

int MyFLA_Set_strict_lower_triangular_to_float( FLA_Obj A, float value );

int MyFLA_Set_strict_lower_triangular_to_double( FLA_Obj A, double value );

int MyFLA_Set_strict_upper_triangular_to_float( FLA_Obj A, float value );

int MyFLA_Set_strict_upper_triangular_to_double( FLA_Obj A, double value );

int MyFLA_Obj_set_to_int( FLA_Obj A, int value );

int MyFLA_Obj_set_to_float( FLA_Obj A, float value );

int MyFLA_Obj_set_to_double( FLA_Obj A, double value );

int MyFLA_Nrm1( FLA_Obj A, FLA_Obj nrm );

int MyFLA_Frob_norm( FLA_Obj A, FLA_Obj nrm );

int MyFLA_Copy_triu( FLA_Obj A, FLA_Obj B );

int MyFLA_Copy( FLA_Obj A, FLA_Obj B );

int MyFLA_Abs( FLA_Obj A );

int MyFLA_Symmetrize_from_lower_matrix( FLA_Obj A );

int MyFLASH_Obj_create( FLA_Datatype datatype, int m, int n, FLA_Obj * A );

int NoFLA_Copy_matrix_d( int m, int n, 
        double * buff_A, int ldim_A, double * buff_B, int ldim_B );

int MyFLA_Generate_random_matrix( FLA_Obj A );

int MyFLA_Generate_spd_matrix( FLA_Obj A );

int MyFLA_Generate_int_lower_triangular( int num, FLA_Obj A );

int MyFLA_Generate_int_upper_triangular( int num, FLA_Obj A );

int MyFLA_Generate_int_matrix( int num, FLA_Obj A );

int MyFLA_Scale_matrix_down( FLA_Obj A );

int MyFLA_Set_matrix_main_diag_from_vector( FLA_Obj v, FLA_Obj A );

int MyFLA_Set_vector_from_matrix_main_diag( FLA_Obj A, FLA_Obj v );

