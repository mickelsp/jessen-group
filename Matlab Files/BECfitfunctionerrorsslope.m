function Z=BECfitfunctionerrorsslope(coeffs,data)
amplitude=coeffs(1);
x0=coeffs(2);
y0=coeffs(3);
xoffset=coeffs(4);
yoffset=coeffs(5);
offset=coeffs(6);
slopex=coeffs(7);
slopey=coeffs(8);

x=data(:,1); %Split the data matrix into x and y vectors
y=data(:,2);
w=data(:,3); %Matrix of errors
numberofpoints=data(:,4);%matrix in which each entry is the number of points 
sqrtnumbpoints=sqrt(numberofpoints(1,1)); %sqrt of number of points

%amplitude=0.36;
%x0=17.91;
%y0=18.99;
%xoffset=91.66;
%yoffset=97.94;
 %offset=0;
 %slopex=0;
%slopey=0;
rotangle=0*pi/180; % It was 30 degree in 84Sr data fitting

         xm=x-xoffset;
         ym=y-yoffset;
         xpsqrd=((cos(rotangle).*xm).^2+(2.*cos(rotangle).*sin(rotangle).*xm.*ym)+(sin(rotangle).*ym).^2);
         ypsqrd=((sin(rotangle).*xm).^2-(2.*cos(rotangle).*sin(rotangle).*xm.*ym)+(cos(rotangle).*ym).^2);


%  Z = (offset+slopex.*(xm)+slopey.*(ym)+...
%      heaviside(1-((xm)./x0).^2-((ym)./y0).^2).*...
%      amplitude.*(1-((xm)./x0).^2-((ym)./y0).^2).^(3/2))./...
%     (w.*sqrtnumbpoints);
 Z = (offset+slopex.*(xm)+slopey.*(ym)+...
     heaviside(1-(xpsqrd./(x0.^2))-(ypsqrd./(y0.^2))).*...
     amplitude.*(1-(xpsqrd./(x0.^2))-(ypsqrd./(y0.^2))).^(3/2))./...
    (w.*sqrtnumbpoints);


%  Z = (offset+slopex.*(xm)+slopey.*(ym)+...
%      heaviside(1-(xpsqrd./x0.^2)-(ypsqrd./y0.^2)).*...
%      amplitude.*(1-((xpsqrd)./x0.^2)-((ypsqrd)./y0.^2)).^(3/2))./...
%      (w.*sqrtnumbpoints); % We're trying to fit z = f(x, y) so compute
%      f(x, y); 


% Z = heaviside(1-((x-xoffset)./x0).^2-((y-yoffset)./y0).^2).*...
%     ((offset+slopex.*(x-xoffset)+slopey.*(y-yoffset)+amplitude.*(1-((x-xoffset)./x0).^2-((y-yoffset)./y0).^2).^(3/2))./...
%     (w.*sqrtnumbpoints)); % We're trying to fit z = f(x, y) so compute f(x, y); 
%     %Heaviside function zeros out everything when Z function becomes complex valued