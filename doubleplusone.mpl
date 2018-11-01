macro(Multiply=LinearAlgebra[Multiply]):
macro(Matrix=LinearAlgebra[Matrix]):
macro(Copy=LinearAlgebra[Copy]):
macro(mInverse=LinearAlgebra[Modular][Inverse]):



# DoublePlusOneLift : Computes a sparse inverse expansion of the matrix A by using Double Plus One Lifting
#
# Input :
#		A : Integer matrix n x n
#		X : Integer which must be relatively prime with the determinant of A
#		n : Dimension of A
#		k : Determines the number of the times to lift and the precision of the inverse, i.e. X^(2^(k+1)-1)
#
# Output:
#		A_0 : It is the inverse of A modulo X (Initialisation of the inverse)
#		R : A table from 0 to k needed for the sparse inverse expansion
#		M : A table from 0 to k - 1 needed for the sparse inverse expansion
#
# All together it will be : (...(A_0(I + R_0 * X_0) + M_0 * X_0^2) * (I + R_1 * X_1) + M_1 * X_1^2)...) = A^(-1) mod X^(2^(k+1)-1
#
DoublePlusOneLift := proc(A, X, n, k)
	A_0 := mInverse(X, A):
	R[0] := map(iquo, 1 - A . A_0, X):
	for i from 0 to (k - 1) do
		R_bar := Multiply(R[i], R[i]):
		M[i] := map(modp, A_0 . R_bar, X):
		R[i+1] := (1 / X) * (R_bar - Multiply(A, M[i]):
	od:
	return A_0, R, M:
end proc:


# XadicRepresentationMatrixCreate : Computes the X-adic representation of a matrix and stores it compactly in a new larger one
#
# Input :
#		A : Input integer matrix n x m
#		X : Integer for the X-adic representation
#		n, m : Dimensions of A
#		p : Length of the X-adic representation to output
#
# Output :
#		A_xadic : Integer matrix n x (m*p) where every vertical slice contains one term from the X-adic representation of A
#
XadicRepresentationMatrixCreate := proc(A, X, n, m, p)
	A_xadic := Matrix(n, m*p):
	A_temp := Copy(A):
	for k to p do
		for i to n do
			for j to m do
				A_xadic[i, j + m * (k - 1)] := modp(A_temp[i, j], X):
				A_temp[i, j] := iquo(A_temp[i, j], X):
			od:
		od:
	od:
	return A_xadic:
end proc:


# XadicRepresentationMatrixCollapse : Computes the actual matrix from its X-adic representation
#
# Input :
#		A_xadic : Input integer matrix n x (m*p) containing the X-adic representation of A
#		X : Integer for the X-adic representation
#		n, m : Dimensions of A to output
#		p : Length of the input X-adic representation
#
# Output :
#		A : Integer matrix n x m collapsed from its X-adic representation
#
XadicRepresentationMatrixCollapse := proc(A_xadic, X, n, m, p)
	A := Matrix(n, m):
	for i to n do
		for j to m do
			sum := 0:
			for k to p do
				sum := sum + A_xadic[i, j + m * (k - 1)] * X^(k - 1)
			od:
			A[i, j] := sum:
		od:
	od:
	return A:
end_proc:


# SwiftRightXadic : Multiplies an X-adic representation by X^swift, i.e., 
#	slides the vertical slices of the matrix swift times at the right
#
# Input :
#		A : Input integer matrix n x (m*p) containing an X-adic representation
#		X : Integer for the X-adic representation
#		n, m : Dimensions of each slice in A
#		p : Length of the X-adic representation, i.e., number of slices
#
# Output :
#		A : Integer matrix n x (m*p) containing the input X-adic representation multiplied by X^swift
#
SwiftRightXadic := proc(A, swift, n, m, p)
	for k from p by -1 to 1 do
		for i to n do
			for j to m do
				# Slide, and set to 0 the first swift slices
				A[i, j + m * (k - 1)] := `if`(k > swift, A[i, j + m * (k - 1 - swift)], 0):
			od:
		od:
	od:
	return A:
end proc:


# CleanUpXadic : Cleans up an X-adic representation, i.e., 
#	scans the matrix, looking for overflowed elements which will reduce mod X, and it will forward the carry to the next slice
#
# Input :
#		A : Input integer matrix n x (m*p) containing an X-adic representation to be cleaned up
#		X : Integer for the X-adic representation
#		n, m : Dimensions of each slice in A
#		p : Length of the X-adic representation, i.e., number of slices
#
# Output :
#		A : Integer matrix n x (m*p) containing the input X-adic representation that has been cleaned up
#
CleanUpXadic := proc(A, X, n, m, p)
	for k to p-1 do
		for i to n do
			for j to m do
				entry := A[i, j + m * (k - 1)]:
				if entry >= X then
					A[i, j + m * (k - 1)] := modp(entry, X):
					A[i, j + m * k] := A[i, j + m * k] + iquo(entry, X):
				fi:
			od:
		od:
	od:
	# Take care of the last slice as well
	for i to n do
		for j to m do
			A[i, j + m * (p - 1)] := modp(A[i, j + m * (p - 1)], X):
		od:
	od:
	return A
end proc:


# ApplyDPOL : Applies the Double Plus One Lifting formula to a matrix B
#
# Input :
#		A_0, R, M : The output from the DoublePlusOneLift procedure, i.e., matrices that form the expansion of the inverse of a matrix A
#		B : Integer matrix n x (m*p) containing an X-adic representation on which the Double Plus One Lifting formula will be applied
#		X : Integer for the X-adic representation
#		n: Dimension of A_0, R[i], M[i], and the number of rows in B
#		m : Number of collumns of each slice in B
#		p : Length of the X-adic representation in B, i.e., number of slices
#		k : Length of the lists R and M
#
# Output :
#		Expansion : Integer matrix n x (m*p) containg the X-adic representation of A^(-1) * B, where A is represented by its sparse inverese expansion here
#
ApplyDPOL := proc(A_0, R, M, B, X, n, m, p, k)
	Expansion := Matrix(n, m*p):
	Factor := B:
	for i from k-1 by -1 to 0 do # k or k-1? See DPOL algorithm
		Expansion := CleanUpXadic(Expansion + Multiply(M[i], SwiftRight(Factor, 2))):
		Factor := CleanUpXadic(Factor + Multiply(R[i], SwiftRight(Factor, 1))):
	Expansion := CleanUpXadic(Expansion + Multiply(A_0, Factor)):
	return Expansion:
#	if k < 0 then
#		return A_0 . B:
#	else
#		return CleanUpXadic(Multiply(M[k], SwiftRight(B, 2)) + 
#			ApplyDPOL(A_0, R, M, CleanUpXadic(B + Multiply(R[k], SwiftRight(B, 1)), X, n, m, p), X, n, m, p, k - 1), X, n, m, p):
#	fi:
end proc:



SolveLinearSystem := proc(A, B, n, m)

	xx := 31:
	k := 100:
	p := 2^k:
	X := [seq(xx^(2^(i+1)-1), i=0..k)]:

	A_0, R, M := DoublePlusOneLift(A, xx, n, k):

	B_xadic := XadicRepresentationMatrixCreate(B, X, n, m, p):

	Solution_xadic := ApplyDPOL(A_0, R, M, B_xadic, X, n, m, p, k):

	return XadicRepresentationMatrixCollapse(Solution_xadic, m, p):

end proc: