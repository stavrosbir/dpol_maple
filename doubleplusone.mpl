macro(RandomMatrix=LinearAlgebra[RandomMatrix]):
macro(HermiteForm=LinearAlgebra[HermiteForm]):
macro(Determinant=LinearAlgebra[Determinant]):
macro(MatrixInverse=LinearAlgebra[MatrixInverse]):
macro(Multiply=LinearAlgebra[Multiply]):
macro(Matrix=LinearAlgebra[Matrix]):
macro(Copy=LinearAlgebra[Copy]):
macro(mMod=LinearAlgebra[Modular][Mod]):
macro(mInverse=LinearAlgebra[Modular][Inverse]):
macro(mCreate=LinearAlgebra[Modular][Create]):
macro(mMultiply=LinearAlgebra[Modular][Multiply]):



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
	A_0 := Inverse(X, A):
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


#
SwiftRightXadic := proc(A, swift, n, m, p)
	for k from p by -1 to swift + 1 do
		for i to n do
			for j to m do
				A[i, j + m * (k - 1)] := A[i, j + m * (k - 1 - swift)]
			od:
		od:
	od:
	return A:
end proc:


#
CleanUpXadic := proc(A, X, n, m, p)
	for k to p do
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
	return A
end proc:


#
ApplyDPOL := proc(A_0, R, M, B, X, n, m, n, p, i)
	if i < 0 then
		return A_0 . B:
	else
		return CleanUpXadic(Multiply(M[i], SwiftRight(B, 2)) + 
			ApplyDPOL(A_0, R, M, CleanUpXadic(B + Multiply(R[i], SwiftRight(B, 1)), X, n, m, p), X, n, m, p, i - 1), X, n, m, p):
	fi:
end proc:



SolveLinearSystem := proc(A, B, n, m, p)

	xx := 31:
	k := 100:
	X := [seq(xx^(2^(i+1)-1), i=0..k)]:

	A_0, R, M := DoublePlusOneLift(A, xx, n, k):

	B_xadic := XadicRepresentationMatrixCreate(B, X, n, m, p):

	Solution_xadic := ApplyDPOL(A_0, R, M, B_xadic, X, n, m, p, k):

	return XadicRepresentationMatrixCollapse(Solution_xadic, m, p):

end proc: