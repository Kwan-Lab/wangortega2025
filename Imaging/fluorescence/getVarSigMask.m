function sigMask =  getVarSigMask(regData) 

if isempty(regData.notSigMask)
            notSig = zeros(size(regData.sigMask,1),size(regData.sigMask,2));
        else
            notSig = combineMasks(regData.notSigMask);
end
     sigMask = ~notSig;
     