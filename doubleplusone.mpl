macro(Multiply=LinearAlgebra[Multiply]):
macro(Copy=LinearAlgebra[Copy]):
macro(mInverse=LinearAlgebra[Modular][Inverse]):



# DoublePlusOneLift : Computes a sparse inverse expansion of the matrix A by using Double Plus One Lifting
#
# Input :
#		A : Integer matrix n x n
#		X : Integer which must be relatively prime with the determinant of A and >= max(10000, 3.61 * n^2 * ||A||)
#		n : Dimension of A
#		k : Determines the number of the times to lift and the precision of the inverse, i.e. X^(2^(k+1)-1)
#
# Precondition : gcd(X, det(A)) = 1 and X >= max(10000, 3.61 * n^2 * ||A||)
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
	R[0] := map(iquo, 1 - Multiply(A, A_0), X):
	for i from 0 to (k - 1) do
		R_bar := Multiply(R[i], R[i]):
		M[i] := map(modp , Multiply(A_0, R_bar), X):
		R[i+1] := map(iquo, R_bar - Multiply(A, M[i]), X):
	od:
	return A_0, eval(R), eval(M):
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
	for i to n do
		for j to m do
			temp := convert(A[i, j], 'base', X):
			for k to min(p, nops(temp)) do
				A_xadic[i, j + m * (k - 1)] := temp[k]:
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
end proc:


# SwiftRightXadic : Multiplies an X-adic representation by X^swift, i.e., 
#	slides the vertical slices of the matrix swift times at the right
#
# Input :
#		A : Input integer matrix n x (m*p) containing an X-adic representation
#		swift : Exponent of X, i.e., A will be multiplied by X^swift
#		n, m : Dimensions of each slice in A
#		p : Length of the X-adic representation, i.e., number of slices
#
# Output :
#		A_slided : Integer matrix n x (m*p) containing the input X-adic representation multiplied by X^swift
#
SwiftRightXadic := proc(A, swift, n, m, p)
	A_slided := Matrix(n, m*p):
	if swift < p then
		A_slided[1..n, swift*m+1..m*p] := A[1..n, 1..m*(p-swift)]:
	fi:
	return A_slided:
end proc:


# SwiftRightAndMultiply : Computes the product of 2 matrices A, B, multiplied by X^swift
#
# Input :
#		A : Integer matrix n x n
#		B : Integer matrix n x (m*p) containing an X-adic representation
#		swift : Exponent of X, i.e., A*B will be multiplied by X^swift
#		n: Dimension of A, and the number of rows in B
#		m : Number of collumns of each slice in B
#		p : Length of the X-adic representation in B, i.e., number of slices
#
# Output :
#		C : Integer matrix n x (m*p) containing the X-adic representation of A * B * X^swift
#
SwiftRightAndMultiply := proc(A, B, swift, n, m, p)
	C := Matrix(n, m*p):
	if swift < p then
		C[1..n, swift*m+1..m*p] := Multiply(A, B[1..n, 1..m*(p-swift)]):
	fi:
	return C:
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
#		A[1..n, 1+m*k..m*(k+1)] := A[1..n, 1+m*k..m*(k+1)] + map(iquo, A[1..n, 1+m*(k-1)..m*k], X):
#		A[1..n, 1+m*(k-1)..m*k] := map(modp, A[1..n, 1+m*(k-1)..m*k], X):
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
#	A[1..n, 1+m*(p-1)..m*p] := map(modp, A[1..n, 1+m*(p-1)..m*p], X):
	for i to n do
		for j to m do
			A[i, j + m * (p - 1)] := modp(A[i, j + m * (p - 1)], X):
		od:
	od:
	return A:
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
# Precondition : p <= 2^(k+1) - 1
#
# Output :
#		Expansion : Integer matrix n x (m*p) containg the X-adic representation of A^(-1) * B, where A is represented by its sparse inverese expansion here
#
ApplyDPOL := proc(A_0, R, M, B, X, n, m, p, k)
	Expansion := Matrix(n, m*p):
	Factor := Copy(B):
	for i from k-1 by -1 to 0 do
		Expansion := CleanUpXadic(Expansion + SwiftRightAndMultiply(M[i], Factor, 2^(i + 2) - 2, n, m, p), X, n, m, p):
		Factor := CleanUpXadic(Factor + SwiftRightAndMultiply(R[i], Factor, 2^(i + 1) - 1, n, m, p), X, n, m, p):
	od:
	Expansion := CleanUpXadic(Expansion + Multiply(A_0, Factor), X, n, m, p):
	return Expansion:
end proc:



# SolveLinearSystem : Solves a linear system A*x=B <=> x=A^(-1)*B over the integers
#
# Input :
#		A : Integer matrix n x n
#		B : Integer matrix n x m
#		n: Dimension of A, and the number of rows in B
#		m : Number of collumns in B
#
# Output :
#		Solution : Integer matrix n x m containg the solution x of the linear system A*x=B
#
SolveLinearSystem := proc(A, B, n, m)

	X := 10007:
	k := 5:
	p := 2^(k+1)-1:

	A_0, R, M := DoublePlusOneLift(A, X, n, k):

	B_xadic := XadicRepresentationMatrixCreate(B, X, n, m, p):

	Solution_xadic := ApplyDPOL(A_0, R, M, B_xadic, X, n, m, p, k):

	Solution := XadicRepresentationMatrixCollapse(Solution_xadic, X, n, m, p):

	return Solution:

end proc:
