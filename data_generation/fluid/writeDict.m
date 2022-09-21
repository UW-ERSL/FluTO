
function writeDict(a)
titleStr1 = "'fill':{'type':'C', 'order':5,\\\n";
formatSpec0 = "'00':{'5':%e,'4':%e,'3':%e,'2':%e,'1':%e,'0':%e},\\\n";
formatSpec1 = "'11':{'5':%e,'4':%e,'3':%e,'2':%e,'1':%e,'0':%e},\\\n";

titleStr2 = "'fill':{'type':'Area', 'order':5,\\\n";
formatSpec2 = "'solid':{'5':%e,'4':%e,'3':%e,'2':%e,'1':%e,'0':%e}";




fileID = fopen('coeffs.txt','w');
fprintf(fileID, titleStr1);
fprintf(fileID, formatSpec0, a(1,:));
fprintf(fileID, formatSpec1, a(2,:));
fprintf(fileID, titleStr2);
fprintf(fileID, formatSpec2, a(3,:));
fclose(fileID);
