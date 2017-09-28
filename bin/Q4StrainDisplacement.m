% This fuction has to be modified in line with the introduced feObj
function B= Q4StrainDisplacement(xi,eta,vertex)
    [~,dsp]=Q4ShapeFn(xi,eta);
    jacobian=dsp*vertex;
    inv_jacobian=inv(jacobian);
    r3=zeros(4,8);
    r1=[1 0 0 0;0 0 0 1;0 1 1 0];
    r2=[inv_jacobian zeros(2,2);zeros(2,2) inv_jacobian];
    r3(1,1:4)=dsp(1,:);
    r3(2,1:4)=dsp(2,:);
    r3(3,5:8)=dsp(1,:);
    r3(4,5:8)=dsp(2,:);
    B=r1*r2*r3;
end
