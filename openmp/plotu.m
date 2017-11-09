writerObj=VideoWriter('test.avi');  %// ����һ����Ƶ�ļ������涯��  
writerObj.FrameRate=3;
open(writerObj);                    %// �򿪸���Ƶ�ļ�  
for i = 0:50:5000  
    f1 = figure(1);
    hold on
%     colormap(gray(256));
    colormap(jet);
    ap=['u', num2str(i)];
    T = readtable(ap);
    area=table2array(T);
    imagesc(1:300,1:300,area(1:300,1:300))%eval: Execute MATLAB expression in text
    axis([1 300 1 300])
    colorbar()
    if(i==0)
        title('t=0.02ms');
    else
        title(['t=' num2str(i*10) 'ms']);
    end
    pause(0.1)
    frame = getframe(f1);            %// ��ͼ�������Ƶ�ļ���  
    writeVideo(writerObj,frame); %// ��֡д����Ƶ  
end
close(writerObj); %// �ر���Ƶ�ļ����  