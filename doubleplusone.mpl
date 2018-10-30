macro(RandomMatrix=LinearAlgebra[RandomMatrix]):
macro(HermiteForm=LinearAlgebra[HermiteForm]):
macro(Determinant=LinearAlgebra[Determinant]):
macro(MatrixInverse=LinearAlgebra[MatrixInverse]):
macro(Multiply=LinearAlgebra[Multiply]):
macro(Matrix=LinearAlgebra[Matrix]):
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
#		B_0 : It is the inverse of A modulo X (Initialisation of the inverse)
#		R : A table from 0 to k needed for the sparse inverse expansion
#		M : A table from 0 to k-1 needed for the sparse inverse expansion
#
# All together it will be : (...(B_0(I + R_0 * X_0) + M_0 * X_0^2) * (I + R_1 * X_1) + M_1 * X_1^2)...) = A^(-1) mod X^(2^(k+1)-1
#
DoublePlusOneLift := proc(A, X, n, k)
	B_0 := Inverse(X, A):
	R[0] := map(iquo, 1 - A . B_0, X):
	for i from 0 to (k-1) do
		R_bar := Multiply(R[i], R[i]):
		M[i] := map(modp, B_0 . R_bar, X):
		R[i+1] := (1 / X) * (R_bar - Multiply(A, M[i]):
	od:
	return B_0, R, M:
end proc:


# XadicRepresentationMatrixCreate : Computes the X-adic representation of a matrix and stores it compactly in a new larger one
#
# Input :
#		A : Input integer matrix n x m
#		X : Integer for the X-adic representation
#		n, m : Dimensions of B
#		p : Length of the X-adic representation to output
#
# Output :
# 		write here
# 
# XadicRepresentationMatrixCreate := proc(A, X, n, m, p)
# 	A_xadic := 
# end proc:


# 
# XadicRepresentationMatrixCollapse := proc(A_xadic, X, n, m, p)
# 	
# end_proc:


# change name
# use cleanup after each mult
# 
ApplyDPOL := proc(B_0, R, M, X, B, i)
	if i < 0 then
		return B_0 . B:
	else
		return M[i] . SwiftRight(B, 2) + ApplyDPOL(B_0, R, M, X, B + R[i] . SwiftRight(B, 1), i - 1):
	end if:
end proc:


SolveLinearSystem := proc(A, B, n, m)
	xx := 31:
	k := 100:
	X := [seq(xx^(2^(i+1)-1), i=0..k)]:
	B_0, R, M := DoublePlusOneLift(A, xx, n, k):
	B_xadic := XadicRepresentationMatrixCreate(B, X, n, m, p):
	Solution_xadic := ApplyDPOL(B_0, R, M, X, B_xadic, k):
	return XadicRepresentationMatrixCollapse(Solution_xadic, m, p):
end proc: