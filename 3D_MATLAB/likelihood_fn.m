function likelihood = likelihood_fn(ObjValue)

if ObjValue > 0
    likelihood = exp(-(ObjValue - (-0) ));
else
    likelihood = exp(-(ObjValue - (-0.5) ));

end