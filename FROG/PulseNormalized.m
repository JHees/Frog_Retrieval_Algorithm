function [P,P_sp] = PulseNormalized(T,P)
% 对于脉冲进行自动归一化，去除零阶相位、一阶相位
    P_sp = fftshift(fft(P));
    [~,region]=PulseMainWidth(F,abs(P_sp),0.5);
    if(region(1)>region(2))
        P = P.*exp(2i*pi*F*T(end));
        P_sp = fftshift(fft(P));
    end
    t_center = sum(T.*abs(P).^2)/sum(abs(P).^2);
    P_sp_offset = P_sp.*exp(2i*pi*F*t_center);
    P_offset = ifft(ifftshift(P_sp_offset));
    f_center = sum(F.*abs(P_sp_offset).^2)/sum(abs(P_sp_offset).^2);
    P = P_offset .*exp(-2i*pi*T*f_center);

    [~, max_ind] = max(abs(P).^2);
    P = P.*exp(-1i*angle(P(max_ind)));
    P_sp = fftshift(fft(P));
end