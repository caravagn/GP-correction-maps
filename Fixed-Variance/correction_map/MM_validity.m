function val = MM_validity(m)
	S_0 = 60;
	kinetics = [2 1 1.5];
	val = m / (S_0+ (kinetics(2) + kinetics(3))/kinetics(1));
return