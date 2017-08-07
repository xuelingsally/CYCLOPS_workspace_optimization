function likelihood = likelihood_fn(ObjValue)

if ObjValue > 0.5
    likelihood = 0;
else
    likelihood = exp(-(ObjValue - (-0.2) ));

end