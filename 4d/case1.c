// ============================================================================
// CONTRIBUTORS
//
// Edoardo Di Napoli
//   Juelich Supercomputing Centre, Institute for Advanced Simulation
//   Forschungszentrum Juelich
//   Wilhelm-Johnen strasse, 52425 Juelich, Germany
//   e.di.napoli@fz-juelich.de
// Diego Fabregat, Paolo Bientinesi
//   AICES, RWTH-Aachen University, 52056-Aachen, Germany
//   {fabregat,pauldj}@aices.rwth-aachen.de
// G. Quintana-Orti
//   Depto. de Ingenieria y Ciencia de Computadores 
//   Universitat Jaume I, 12.080 Castellon, Spain
//   gquintan@icc.uji.es
//
// ============================================================================
#include <math.h>
#include "FLAME.h"
#include "MyFLA_Utils.h"
#include "case1.h"


// ============================================================================
void compute_case1ref( int size_a, int size_b, int size_c, int size_d,
         int size_i, int size_j, FLA_Obj cb_A, FLA_Obj cb_B, FLA_Obj cb_C, 
         int print_data ) {
  int      size_ab, size_abc, size_ia, size_iaj, size_ic, size_icj,
           iter_a, iter_b, iter_c, iter_d, iter_i, iter_j;
  double   * buff_cb_A, * buff_cb_B, * buff_cb_C, ci;

  // Some initializations.
  buff_cb_A = ( double * ) FLA_Obj_buffer_at_view( cb_A );
  buff_cb_B = ( double * ) FLA_Obj_buffer_at_view( cb_B );
  buff_cb_C = ( double * ) FLA_Obj_buffer_at_view( cb_C );
  
  size_ab  = size_a * size_b;
  size_abc = size_a * size_b * size_c;

  size_ia  = size_i * size_a;
  size_iaj = size_i * size_a * size_j;

  size_ic  = size_i * size_c;
  size_icj = size_i * size_c * size_j;

  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " cb_A_i = [ ", cb_A, "%le", " ];" );
    FLA_Obj_show( " cb_B_i = [ ", cb_B, "%le", " ];" );
    FLA_Obj_show( " cb_C_i = [ ", cb_C, "%le", " ];" );
  }

  // Perform computation.
  for( iter_a = 0; iter_a < size_a; iter_a++ ) {
    for( iter_b = 0; iter_b < size_b; iter_b++ ) {
      for( iter_c = 0; iter_c < size_c; iter_c++ ) {
        for( iter_d = 0; iter_d < size_d; iter_d++ ) {
          ci = 0.0;
          for( iter_i = 0; iter_i < size_i; iter_i++ ) {
            for( iter_j = 0; iter_j < size_j; iter_j++ ) {
              ci += buff_cb_A[ iter_i + iter_a * size_i + iter_j * size_ia + 
                                   iter_b * size_iaj ] *
                    buff_cb_B[ iter_i + iter_c * size_i + iter_j * size_ic + 
                                   iter_d * size_icj ];
            }
          }
          buff_cb_C[ iter_a + iter_b * size_a + iter_c * size_ab + 
                         iter_d * size_abc ] = ci;
        }
      }
    }
  }

  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " cb_A_f = [ ", cb_A, "%le", " ];" );
    FLA_Obj_show( " cb_B_f = [ ", cb_B, "%le", " ];" );
    FLA_Obj_show( " cb_C_f = [ ", cb_C, "%le", " ];" );
  }
}


// ============================================================================
void compute_case1a( int size_a, int size_b, int size_c, int size_d,
         int size_i, int size_j, FLA_Obj cb_A, FLA_Obj cb_B, FLA_Obj cb_C, 
         int print_data ) {
  int      size_ab, size_abc, size_ia, size_iaj, size_ic, size_icj,
           iter_a, iter_b, iter_c, iter_d, iter_i, iter_j;
  double   * buff_cb_A, * buff_cb_B, * buff_cb_C, ci, * ptr_ai, * ptr_bi;

  // Some initializations.
  buff_cb_A = ( double * ) FLA_Obj_buffer_at_view( cb_A );
  buff_cb_B = ( double * ) FLA_Obj_buffer_at_view( cb_B );
  buff_cb_C = ( double * ) FLA_Obj_buffer_at_view( cb_C );
  
  size_ab  = size_a * size_b;
  size_abc = size_a * size_b * size_c;

  size_ia  = size_i * size_a;
  size_iaj = size_i * size_a * size_j;

  size_ic  = size_i * size_c;
  size_icj = size_i * size_c * size_j;

  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " cb_A_i = [ ", cb_A, "%le", " ];" );
    FLA_Obj_show( " cb_B_i = [ ", cb_B, "%le", " ];" );
    FLA_Obj_show( " cb_C_i = [ ", cb_C, "%le", " ];" );
  }

  // Perform computation.
  for( iter_a = 0; iter_a < size_a; iter_a++ ) {
    for( iter_b = 0; iter_b < size_b; iter_b++ ) {
      for( iter_c = 0; iter_c < size_c; iter_c++ ) {
        for( iter_d = 0; iter_d < size_d; iter_d++ ) {
          ci = 0.0;
          for( iter_j = 0; iter_j < size_j; iter_j++ ) {
            ptr_ai = & buff_cb_A[ 0 + iter_a * size_i + iter_j * size_ia + 
                                      iter_b * size_iaj ];
            ptr_bi = & buff_cb_B[ 0 + iter_c * size_i + iter_j * size_ic + 
                                      iter_d * size_icj ];
            for( iter_i = 0; iter_i < size_i; iter_i++ ) {
              ci += ( * ptr_ai ) * ( * ptr_bi );
              ptr_ai++;
              ptr_bi++;
            }
          }
          buff_cb_C[ iter_a + iter_b * size_a + iter_c * size_ab + 
                         iter_d * size_abc ] = ci;
        }
      }
    }
  }

  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " cb_A_f = [ ", cb_A, "%le", " ];" );
    FLA_Obj_show( " cb_B_f = [ ", cb_B, "%le", " ];" );
    FLA_Obj_show( " cb_C_f = [ ", cb_C, "%le", " ];" );
  }
}


// ============================================================================
void compute_case1b( int size_a, int size_b, int size_c, int size_d,
         int size_i, int size_j, FLA_Obj cb_A, FLA_Obj cb_B, FLA_Obj cb_C, 
         int print_data ) {
  FLA_Obj  slice_A, slice_B, slice_C;
  int      datatype, size_ab, size_abc, size_ia, size_iaj, size_ic, size_icj,
           iter_a, iter_b, iter_c, iter_d, iter_i, iter_j;
  size_t   idx_A, idx_B, idx_C;
  double   * buff_cb_A, * buff_cb_B, * buff_cb_C;

  // Some initializations.
  datatype  = FLA_Obj_datatype( cb_A );
  buff_cb_A = ( double * ) FLA_Obj_buffer_at_view( cb_A );
  buff_cb_B = ( double * ) FLA_Obj_buffer_at_view( cb_B );
  buff_cb_C = ( double * ) FLA_Obj_buffer_at_view( cb_C );
  
  size_ab  = size_a * size_b;
  size_abc = size_a * size_b * size_c;

  size_ia  = size_i * size_a;
  size_iaj = size_i * size_a * size_j;

  size_ic  = size_i * size_c;
  size_icj = size_i * size_c * size_j;

  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " cb_A_i = [ ", cb_A, "%le", " ];" );
    FLA_Obj_show( " cb_B_i = [ ", cb_B, "%le", " ];" );
    FLA_Obj_show( " cb_C_i = [ ", cb_C, "%le", " ];" );
  }

  // Prepare temporal slices without buffer.
  FLA_Obj_create_without_buffer( datatype, size_i, size_a, & slice_A );
  FLA_Obj_create_without_buffer( datatype, size_i, size_c, & slice_B );
  FLA_Obj_create_without_buffer( datatype, size_a, size_c, & slice_C );

  // Perform computation.
  for( iter_b = 0; iter_b < size_b; iter_b++ ) {
   
    for( iter_d = 0; iter_d < size_d; iter_d++ ) {

      iter_a = 0;
      iter_c = 0;
      iter_i = 0;

      idx_C = ( ( size_t ) iter_a ) + 
              ( ( size_t ) iter_b * size_a ) + 
              ( ( size_t ) iter_c * size_ab ) + 
              ( ( size_t ) iter_d * size_abc );
      FLA_Obj_attach_buffer( & buff_cb_C[ idx_C ], 1, size_ab, & slice_C );
      MyFLA_Obj_set_to_zero( slice_C );
 
      for( iter_j = 0; iter_j < size_j; iter_j++ ) {
        idx_A = ( ( size_t ) iter_i ) +
                ( ( size_t ) iter_a * size_i ) +
                ( ( size_t ) iter_j * size_ia ) +
                ( ( size_t ) iter_b * size_iaj );
        FLA_Obj_attach_buffer( & buff_cb_A[ idx_A ], 1, size_i, & slice_A );

        idx_B = ( ( size_t ) iter_i ) +
                ( ( size_t ) iter_c * size_i ) +
                ( ( size_t ) iter_j * size_ic ) +
                ( ( size_t ) iter_d * size_icj );
        FLA_Obj_attach_buffer( & buff_cb_B[ idx_B ], 1, size_i, & slice_B );

        FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, 
                  FLA_ONE, slice_A, slice_B, FLA_ONE, slice_C );
      }
    }
  }

  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " cb_A_f = [ ", cb_A, "%le", " ];" );
    FLA_Obj_show( " cb_B_f = [ ", cb_B, "%le", " ];" );
    FLA_Obj_show( " cb_C_f = [ ", cb_C, "%le", " ];" );
  }

  // Remove temporal slices.
  FLA_Obj_free_without_buffer( & slice_A );
  FLA_Obj_free_without_buffer( & slice_B );
  FLA_Obj_free_without_buffer( & slice_C );
}


