

# with(LinearAlgebra):
# with(LinearAlgebra[Modular]):

DoublePlusOneLift := proc(A, X, n, k)
	B_0 := Inverse(X, A):
	R[0] := map(iquo, 1 - A . B_0, X):
	for i from 0 to (k-1) do
		R_bar := R[i] . R[i]:
		M[i] := map(modp, B_0 . R_bar, X):
		R[i+1] := (1 / X) * (R_bar - A . M[i]):
	od:
	return B_0, R, M:
end proc:

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
	xx := 1009:
	k := 100:
	X := [seq(xx^(2^(i+1)-1), i=0..k)]:
	B_0, R, M := DoublePlusOneLift(A, xx, n, k):
	B_xadic := XadicRepresentationMatrixCreate(B, m, X[k]):
	Solution_xadic := ApplyDPOL(B_0, R, M, X, B_xadic, k):
	return XadicRepresentationMatrixCollapse(Solution_xadic, m, X[k]):
end proc: