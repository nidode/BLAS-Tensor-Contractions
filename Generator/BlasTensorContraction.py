#!/usr/bin/env python

# TO-DO
# - Get tensor name as input (instead of default T1 T2)
# - Transposition (TransA, TransB)
# - Extract layer for actual executable Matlab code

import sys

def assert_non_repeated_indices( indices ):
	occurrences = set()
	for idx in indices:
		if idx in occurrences:
			print "Multiple occurrences of index %s. Aborting..."
			sys.exit()
		occurrences.add( idx )

def contraction_formula( T1, T2, for_all_indices, sum_indices, op_label ):
	#print for_all_indices
	#print sum_indices
	def str_for_all( for_all_indices ):
		if for_all_indices:
			return " ".join(["Forall_%s" % i for i in for_all_indices ]) + " "
		else:
			return ""
	def str_sum( sum_indices ):
		if sum_indices:
			return " ".join(["Sum_%s" % i for i in sum_indices ]) + " "
		else:
			return ""
	def str_sliced_tensor( T_id, T_sliced ):
		if T_sliced:
			return T_id + "_" + "".join(T_sliced)
		else:
			return T_id
	sliced_indices = for_all_indices + sum_indices
	T1_sliced = [ idx for idx in T1.indices if idx in sliced_indices ]
	T2_sliced = [ idx for idx in T2.indices if idx in sliced_indices ]
	return "%s%s%s %s  %s" % ( str_for_all(for_all_indices), \
	                           str_sum(sum_indices), \
					           str_sliced_tensor( "T1", T1_sliced ), \
					           str_sliced_tensor( "T2", T2_sliced ), \
							   op_label )

class Tensor:
	_indent = "  "
	def __init__(self, indices):
		assert( len(indices) > 0 )
		assert_non_repeated_indices( indices )
		self.indices = indices

	def __mul__(self, T2):
		#contracted = set(self.indices).intersection(set(T2.indices))
		contracted = list(set(self.indices).intersection(set(T2.indices)))
		assert( len(contracted) > 0 )
		#T1_free = set(self.indices).difference(contracted)
		T1_free = [idx for idx in self.indices if idx not in contracted]
		#T2_free = set(T2.indices).difference(contracted)
		T2_free = [idx for idx in T2.indices if idx not in contracted]
		#print contracted
		#print T1_free
		#print T2_free


		# Class 1
		if len(T1_free) == 0 and len(T2_free) == 0:
			#print "Class 1"
			for_all_indices = []
			sum_indices = self.indices[1:]
			return contraction_formula( self, T2, for_all_indices, sum_indices, "(DOT)" )
		# Class 2
		elif len(T1_free) > 0 and len(T2_free) == 0 or \
		     len(T1_free) == 0 and len(T2_free) > 0:
			#print "Class 2"
			if len(T1_free) > 0:
				for_all_indices = T1_free[1:]
				if self.indices[0] in contracted:
					sum_indices = [ idx for idx in contracted if idx != self.indices[0] ]
				else:
					sum_indices = contracted[1:]
			else: # T2 the one with free indices
				for_all_indices = T2_free[1:]
				if T2.indices[0] in contracted:
					sum_indices = [ idx for idx in contracted if idx != T2.indices[0] ]
				else:
					sum_indices = contracted[1:]
			return contraction_formula( self, T2, for_all_indices, sum_indices, "(GEMV)" )
		# Class 3
		elif len(T1_free) > 0 and len(T2_free) > 0:
			# Class 3.1
			if self.indices[0] in contracted and \
				T2.indices[0] in contracted and \
				self.indices[0] != T2.indices[0]:
				#print "Class 3.1"
				# Take stride 1 of T1 as part of the gemm
				# and thus copy-transpose slices of T2
				# As free indices for gemm we pick te first of each tensor
				for_all_indices = T1_free[1:] + T2_free[1:]
				sum_indices = [ idx for idx in contracted if idx != self.indices[0] ]
				#
				return contraction_formula( self, T2, for_all_indices, sum_indices, "(COPY+GEMM)" )
			# Class 3.2
			else:
				#print "Class 3.2"
				# Take first index (stride 1) of each tensor
				T1_str1 = self.indices[0]
				T2_str1 = T2.indices[0]
				# If same (contracted index) take one free each for the gemm
				# and slice all the others (free and contracted)
				if T1_str1 == T2_str1:
					for_all_indices = T1_free[1:] + T2_free[1:]
					sum_indices = [ idx for idx in contracted if idx != T1_str1 ]
				else:
					# If both are free, take the first contracted for the gemm
					if T1_str1 in T1_free and T2_str1 in T2_free:
						for_all_indices = T1_free[1:] + T2_free[1:]
						sum_indices = contracted[1:]
					# T1_str1 is contracted, take the first free index of T1 for the gemm
					elif T1_str1 in contracted and T2_str1 in T2_free:
						for_all_indices = T1_free[1:] + T2_free[1:]
						sum_indices = [ idx for idx in contracted if idx != T1_str1 ]
					# T2_str1 is the contracted one, take its first free index for the gemm
					else:
						for_all_indices = T1_free[1:] + T2_free[1:]
						sum_indices = [ idx for idx in contracted if idx != T2_str1 ]
				return contraction_formula( self, T2, for_all_indices, sum_indices, "(GEMM)" )

				

if __name__ == "__main__":
	def check_test( output, expected ):
		#print output
		#print expected
		if output == expected:
			print " [OK]"
		else:
			print " [FAILED]"

	###########
	# Class 1
	###########
	print "Testing CLASS 1:\n"

	#print "** T1_i * T2_i"
	T1 = Tensor(["i"])
	T2 = Tensor(["i"])
	print "Single contraction...",
	check_test( T1 * T2, "T1 T2  (DOT)" )

	#print "** T1_ijk * T2_ijk"
	T1 = Tensor(["i", "j", "k"])
	T2 = Tensor(["i", "j", "k"])
	print "Triple contraction ordered...",
	check_test( T1 * T2, "Sum_j Sum_k T1_jk T2_jk  (DOT)" )

	#print "** T1_ijk * T2_kji"
	T1 = Tensor(["i", "j", "k"])
	T2 = Tensor(["k", "j", "i"])
	print "Triple contraction shuffled...",
	check_test( T1 * T2, "Sum_j Sum_k T1_jk T2_kj  (DOT)" )

	print

	###########
	# Class 2
	###########
	print "Testing CLASS 2:\n"

	T1 = Tensor(["j", "i"])
	T2 = Tensor(["j"])
	print "Both stride-1 contracted and same, T1 has free indices (1c1f/1c0f)...",
	check_test( T1 * T2, "T1 T2  (GEMV)" )

	T1 = Tensor(["j", "i", "l"])
	T2 = Tensor(["j"])
	print "Both stride-1 contracted and same, T1 has free indices (1c2f/1c0f)...",
	check_test( T1 * T2, "Forall_l T1_l T2  (GEMV)" )

	T1 = Tensor(["j", "k", "i", "l"])
	T2 = Tensor(["i", "j"])
	print "Both stride-1 contracted and distinct, T1 has free indices (2c2f/2c0f)...",
	check_test( T1 * T2, "Forall_l Sum_i T1_il T2_i  (GEMV)" )

	T1 = Tensor(["i", "j"])
	T2 = Tensor(["j"])
	print "T1 stride-1 is free, T1 has free indices (1c1f/1c0f)...",
	check_test( T1 * T2, "T1 T2  (GEMV)" )

	T1 = Tensor(["i", "j", "k", "l"])
	T2 = Tensor(["j", "l"])
	print "T1 stride-1 is free, T1 has free indices (2c2f/2c0f)...",
	check_test( T1 * T2, "Forall_k Sum_l T1_kl T2_l  (GEMV)" )

	T1 = Tensor(["j"])
	T2 = Tensor(["j", "i"])
	print "Both stride-1 contracted and same, T2 has free indices (1c0f/1c1f)...",
	check_test( T1 * T2, "T1 T2  (GEMV)" )

	T1 = Tensor(["j"])
	T2 = Tensor(["j", "i", "l"])
	print "Both stride-1 contracted and same, T2 has free indices (1c0f/1c2f)...",
	check_test( T1 * T2, "Forall_l T1 T2_l  (GEMV)" )

	T1 = Tensor(["i", "j"])
	T2 = Tensor(["j", "k", "i", "l"])
	print "Both stride-1 contracted and distinct, T2 has free indices (2c0f/2c2f)...",
	check_test( T1 * T2, "Forall_l Sum_i T1_i T2_il  (GEMV)" )

	T1 = Tensor(["j"])
	T2 = Tensor(["i", "j"])
	print "T2 stride-1 is free, T2 has free indices (1c0f/1c1f)...",
	check_test( T1 * T2, "T1 T2  (GEMV)" )

	T1 = Tensor(["j", "l"])
	T2 = Tensor(["i", "j", "k", "l"])
	print "T2 stride-1 is free, T2 has free indices (2c0f/2c2f)...",
	check_test( T1 * T2, "Forall_k Sum_l T1_l T2_kl  (GEMV)" )

	print

	###########
	# Class 3.1
	###########
	print "Testing CLASS 3.1:\n"

	T1 = Tensor(["i", "j", "l"])
	T2 = Tensor(["l", "i", "m"])
	print "Both stride-1 contracted and distinct (2c1f/2c1f)...",
	check_test( T1 * T2, "Sum_l T1_l T2_l  (COPY+GEMM)" )

	#print "** T1_ijkl * T2_limn"
	T1 = Tensor(["i", "j", "l", "n"])
	T2 = Tensor(["l", "i", "m", "n", "k"])
	print "Both stride-1 contracted and distinct (3c1f/3c2f)...",
	check_test( T1 * T2, "Forall_k Sum_l Sum_n T1_ln T2_lnk  (COPY+GEMM)" )

	#print "** T1_ijl * T2_lki"
	#T1 = Tensor(["i", "j", "l"])
	#T2 = Tensor(["l", "k", "i"])
	#T1 * T2
	#print "Sum_l T1_l T2_l?"
	#print

	print

	###########
	# Class 3.2
	###########
	print "Testing CLASS 3.2:\n"

	T1 = Tensor(["i", "j"])
	T2 = Tensor(["i", "m"])
	print "Both stride-1 contracted and same (1c1f/1c1f)...",
	check_test( T1 * T2, "T1 T2  (GEMM)" )

	T1 = Tensor(["i", "j", "k"])
	T2 = Tensor(["i", "l"])
	print "Both stride-1 contracted and same (1c2f/1c1f)...",
	check_test( T1 * T2, "Forall_k T1_k T2  (GEMM)" )

	T1 = Tensor(["i", "j"])
	T2 = Tensor(["i", "k", "l"])
	print "Both stride-1 contracted and same (1c1f/1c2f)...",
	check_test( T1 * T2, "Forall_l T1 T2_l  (GEMM)" )

	#print "** T1_ijkl * T2_ilmn"
	T1 = Tensor(["i", "j", "k", "l"])
	T2 = Tensor(["i", "l", "m", "n"])
	print "Both stride-1 contracted and same (2c2f/2c2f)...",
	check_test( T1 * T2, "Forall_k Forall_n Sum_l T1_kl T2_ln  (GEMM)" )

	T1 = Tensor(["i", "k"])
	T2 = Tensor(["m", "i"])
	print "T1 stride-1 contracted (1c1f/1c1f)...",
	check_test( T1 * T2, "T1 T2  (GEMM)" )

	T1 = Tensor(["i", "k", "l"])
	T2 = Tensor(["m", "i", "j"])
	print "T1 stride-1 contracted (1c2f/1c2f)...",
	check_test( T1 * T2, "Forall_l Forall_j T1_l T2_j  (GEMM)" )

	#print "** T1_iklj * T2_mlin"
	T1 = Tensor(["i", "k", "l", "j"])
	T2 = Tensor(["m", "l", "i", "n"])
	print "T1 stride-1 contracted (2c2f/2c2f)...",
	check_test( T1 * T2, "Forall_j Forall_n Sum_l T1_lj T2_ln  (GEMM)" )

	T1 = Tensor(["m", "i"])
	T2 = Tensor(["i", "k"])
	print "T2 stride-1 contracted (1c1f/1c1f)...",
	check_test( T1 * T2, "T1 T2  (GEMM)" )

	T1 = Tensor(["m", "i", "j"])
	T2 = Tensor(["i", "k", "l"])
	print "T2 stride-1 contracted (1c2f/1c2f)...",
	check_test( T1 * T2, "Forall_j Forall_l T1_j T2_l  (GEMM)" )

	#print "** T1_kilj * T2_lmin"
	T1 = Tensor(["k", "i", "l", "j"])
	T2 = Tensor(["l", "m", "i", "n"])
	print "T2 stride-1 contracted (2c2f/2c2f)...",
	check_test( T1 * T2, "Forall_j Forall_n Sum_i T1_ij T2_in  (GEMM)" )

	T1 = Tensor(["k", "i"])
	T2 = Tensor(["n", "i"])
	print "Both stride-1 free (1c1f/1c1f)...",
	check_test( T1 * T2, "T1 T2  (GEMM)" )

	T1 = Tensor(["k", "i", "l"])
	T2 = Tensor(["n", "i"])
	print "Both stride-1 free (1c2f/1c1f)...",
	check_test( T1 * T2, "Forall_l T1_l T2  (GEMM)" )

	T1 = Tensor(["k", "i"])
	T2 = Tensor(["n", "i", "l"])
	print "Both stride-1 free (1c1f/1c2f)...",
	check_test( T1 * T2, "Forall_l T1 T2_l  (GEMM)" )

	#print "** T1_kilj * T2_nmil"
	T1 = Tensor(["k", "i", "l", "j", "p"])
	T2 = Tensor(["n", "m", "i", "l", "q"])
	print "Both stride-1 free (2c3f/2c3f)...",
	check_test( T1 * T2, "Forall_j Forall_p Forall_m Forall_q Sum_l T1_ljp T2_mlq  (GEMM)" )
