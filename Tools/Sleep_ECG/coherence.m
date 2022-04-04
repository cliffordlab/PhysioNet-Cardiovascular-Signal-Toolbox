function co = coherence(x,y,window,nfft)
% x must be power of 2, e.g. 1024

seg=length(x)/2;
overlap=seg/2;
for i=1:3
    x1=x(1+(i-1)*overlap:(i-1)*overlap+seg);
    y1=y(1+(i-1)*overlap:(i-1)*overlap+seg);
    cross_sp(i,:)=cpsd(x1,y1,window,[],nfft);
    h=spectrum.welch;
    if isempty(window)
        window=hanning(length(x1));
    end
    x1=x1'.*window;
    y1=y1'.*window;
    Rx1=psd(h,x1,'NFFT',nfft);
    Ry1=psd(h,y1,'NFFT',nfft);
    Rx(i,:)=Rx1.Data;
    Ry(i,:)=Ry1.Data;
end
m_cross_sp=mean(cross_sp);
m_Rx=mean(Rx);
m_Ry=mean(Ry);
co=(abs(m_cross_sp).^2)./(m_Rx.*m_Ry);
end
