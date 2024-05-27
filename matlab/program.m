%ΟΝΟΜΑΤΕΠΩΝΥΜΟ: Αλεξοπούλου Γεωργία
%ΑΜ: 03120164

%ΑΣΚΗΣΗ 1
    
    % Ερώτημα α
    % εφόσον το αρχείο έχει ήδη ηχογραφηθεί, ο κώδικας για το συγκεκριμένο
    % ερώτημα θα παραμείνει ως σχόλιο, προκειμένου το πρόγραμμα να "τρέξει"
    % κανονικά
    %Fs = 8000; Channels = 1 ; bits = 16;
    %r = audiorecorder(Fs, bits, Channels);
    %duration = 1.5;
    %disp('Start speaking');
    %recordblocking(r, duration);
    %disp('Stop speaking'); 
    %x = getaudiodata(r);
    %
    %filename = 'gina.wav';
    %audiowrite(filename, x, Fs);
    
    Fs = 8000; Channels = 1 ; bits = 16;
    rootdirectory = 'C:\Users\Gina\OneDrive\Desktop\ΜΑΘΗΜΑΤΑ\ΣΧΟΛΗ\Σήματα και Συστήματα\Matlab\03120164';
    [x, Fs] = audioread(fullfile(rootdirectory, 'gina.wav'));
    
    % Ερώτημα β
    t1 = 0 : 1/Fs : (length(x)-1)/Fs;
    figure(1)
    subplot(2,1,1), plot(t1, x, 'LineWidth', 1.5);
    xlabel('time (sec)'); ylabel('Amplitude');
    title('Time Domain Plot of the Recorded Signal');
    sound(x,Fs);
    pause(2)
    
    %επιλέγεται το τμήμα μεταξύ των 0.95-1 sec
    y= x((951*10^(-3))*Fs:(1000*10^(-3))*Fs);
    t2 = 0 : 1/Fs : (length(y)-1)/Fs;
    subplot(2,1,2), plot(t2, y, 'LineWidth', 1.5);
    xlabel('time (sec)'); ylabel('Amplitude');
    title('50msec frame window');
    sound(y,Fs);
    pause(2)
    
    % Ερώτημα γ
    M = 200;
    n= 0:200;
    w = 0.54- 0.46*cos(2*pi*n/M);
    
    N = normalize(x)/5;
    
    figure(2)
    subplot(3,1,1), plot(t1,N);
    title('N');
    
    E = conv(N.^2,w,"valid");
    subplot(3,1,2), plot(E);
    title('E');
    
    subplot(3,1,3), plot(N); hold on; plot(E./max(E), LineWidth = 1); title('w[n]'); hold off;
    
    % Ερώτημα δ
    four=fft(y,1024)./length(y);
    figure(3)
    subplot(2,1,1)
    plot(abs(four))
    title('Time-Domain Signal')
    
    subplot(2,1,2)
    plot(20*log(abs(four(2:end))))
    ylabel('Magnitude (In db)')

% ΆΣΚΗΣΗ 2.1
    % Αρχικά υπολογίζω τα ερωτήματα α-γ για Ν = 2
    % Έρώτημα α
    N1 = 2;
    s1 = N1 + 1;
    a1 = [1 0 0];
    b1 = (1/s1)*ones(1, s1);
    
    % Ερώτημα β
    figure(4)
    freqz(b1,a1)
    
    % Ερώτημα γ
    [z, p, k] = tf2zp(b1,a1);
    figure(5)
    zplane(z, p)
    
    % Υπολογίζω τα ερωτήματα α-γ για Ν = 5
    % Έρώτημα α
    N2 = 5;
    s2 = N2 + 1;
    a2 = [1 0 0 0 0 0];
    b2 = (1/s2)*ones(1, s2);
    
    % Ερώτημα β
    figure(6)
    freqz(b2,a2)
    
    % Ερώτημα γ
    [z, p, k] = tf2zp(b2,a2);
    figure(7)
    zplane(z, p)
    
    % Υπολογίζω τα ερωτήματα α-γ για Ν = 10
    % Έρώτημα α
    N3 = 10;
    s3 = N3 + 1;
    a3 = [1 0 0 0 0 0 0 0 0 0 0];
    b3 = (1/s3)*ones(1, s3);
    
    % Ερώτημα β
    figure(8)
    freqz(b3,a3)
    figure(9)
    
    % Ερώτημα γ
    [z, p, k] = tf2zp(b3,a3);
    zplane(z, p)
    
    % Συνεχίζουμε στο ερώτημα δ
    Wc = 0.1;
    % Για n = 2
    [B1, A1] = butter(2, Wc);
    figure(10)
    subplot(2,1,1)
    freqz(B1,A1)
    subplot(2,1,2)
    [z, p, k] = tf2zp(B1,A1);
    zplane(z, p)
    % Για n = 8
    [B2, A2] = butter(8, Wc);
    figure(11)
    subplot(2,1,1)
    freqz(B2,A2)
    subplot(2,1,2)
    [z, p, k] = tf2zp(B2,A2);
    zplane(z, p)

% Άσκηση 2.2
    % Ερώτημα α
    [x1, Fs] = audioread("cello_note.wav");
    sound(x1, Fs)
    pause(2)
    figure(12)
    subplot(3,1,1)
    plot(x1)
    y = fft(x1)
    subplot(3,1,2)
    plot(y)
    subplot(3,1,3)
    plot(abs(y))
    
    %Ερώτημα β
    [x2, Fs] = audioread("cello_note_noisy.wav");
    sound(x2, Fs)
    pause(2)
    figure(13)
    subplot(3,1,1)
    plot(x2)
    subplot(3,1,2)
    y = fft(x2)
    plot(y)
    subplot(3,1,3)
    plot(abs(y))
    
    %Ερώτημα γ
    %N = 2
    [b, a] = butter(2,0.16);
    denoisation = filter(b, a, x2);
    fft_denoisation = fft(denoisation, 44100);
    figure(34)
    subplot(2,1,1), plot(abs(fft_denoisation))
    title('Butterworth 2ης τάξης - Ν = 2')
    
    %Ν = 12
    [b, a] = butter(12,0.16);
    denoisation_again = filter(b, a, x2);
    fft_denoisation_again = fft(denoisation_again, 44100);
    subplot(2,1,2), plot(abs(fft_denoisation_again))
    title('Butterworth 12ης τάξης - Ν = 12')

    %Σύγκριση ήχων
    sound(denoisation);
    pause(4)
    sound(denoisation_again);
    pause(2)

    %Ερώτημα δ
    N = 10;
    s = N + 1;
    
    a = [1 0 0 0 0 0 0 0 0 0 0];
    b = (1/s)*ones(1,s);
    
    audio_filtered = filter(b, a, x);
    
    pause(1);
    sound(x2, Fs);
    pause(2);
    sound(audio_filtered, Fs);
    pause(2);
    sound(x1, Fs);
    
    figure(14)
    subplot(3,1,1)
    plot(x);
    subplot(3,1,2)
    fft4 = fft(audio_filtered);
    plot(fft4);
    subplot(3,1,3)
    plot(abs(fft4));

% ΆΣΚΗΣΗ 3
    % Ερώτημα α
    % Η συνάρτηση έχει ορισθεί στο αρχείο "resonator.m"
    
    % Ερώτημα β
    fr = 200;
    fs = 400; 
    d = zeros(1, 400);
    d(1) = 1; 
    % Θεωρούμε Ω = π
    % Για r = 0.95
    r = 0.95;
    IR = resonator(d, fr, r, fs);
    IRF = fft(IR);
    figure(15);
    subplot(2, 1, 1);
    time = linspace(0, 1, 400);
    plot(time, IR);
    title('r = 0.95');
    subplot(2, 1, 2);
    freq = linspace(0, 400, 400);
    plot(freq, abs(IRF));
    
    % Για r = 0.7
    r = 0.7;
    IR = resonator(d, fr, r, fs);
    IRF = fft(IR);
    figure(16);
    subplot(2, 1, 1);
    time = linspace(0, 1, 400);
    plot(time, IR);
    title('r = 0.7');
    subplot(2, 1, 2);
    freq = linspace(0, 400, 400);
    plot(freq, abs(IRF));
    
    % Για r = 1
    r = 1;
    IR = resonator(d, fr, r, fs);
    IRF = fft(IR);
    figure(17);
    subplot(2, 1, 1);
    time = linspace(0, 1, 400);
    plot(time, IR);
    title('r = 1');
    subplot(2, 1, 2);
    freq = linspace(0, 400, 400);
    plot(freq, abs(IRF));
    
    % Για r = 1.2
    r = 1.2;
    IR = resonator(d, fr, r, fs);
    IRF = fft(IR);
    figure(18);
    subplot(2, 1, 1);
    time = linspace(0, 1, 400);
    plot(time, IR);
    title('r = 1.2');
    subplot(2, 1, 2);
    freq = linspace(0, 400, 400);
    plot(freq, abs(IRF));
    
    % Ερώτημα γ
    fs = 8000; 
    fr1 = 500;
    fr2 = 1500;
    fr3 = 2500;
    d = zeros(1, 8000);
    d(1) = 1;
    r = 0.95;
    IR1 = resonator(d, fr1, r, fs);
    IR2 = resonator(IR1, fr2, r, fs);
    IR3 = resonator(IR2, fr3, r, fs);
    TotalFR = fft(IR3);
    freq = linspace(0, 8000, 8000);
    figure(19);
    freqz(freq, abs(TotalFR)); 
    
    % Ερώτημα δ
    r = 0.95;
    fs = 8000; 
    fr1 = 500;
    fr2 = 1500;
    fr3 = 2500;
    ImpulseTrain = zeros(1, 1601); 
    ImpulseTrain(1:80:end) = 1; 
    IT1 = resonator(ImpulseTrain, fr1, r, fs);
    IT2 = resonator(IT1, fr2, r, fs);
    IT3 = resonator(IT2, fr3, r, fs);
    TrainRF= fft(IT3);
    % Μέτρο εξόδου DFT 
    figure(20);
    freq = linspace(0, 8000, 1601);
    plot(freq, abs(TrainRF));
    num = [1, -1];
    denum = [1, 0];
    FirstDif = filter(num, denum, IT3); %υπολογισμός διαφοράς x[n]-x[n-1]
    % Σχεδιασμός x[n]
    figure(21);
    subplot(2, 1, 1);
    time = linspace(0, 0.2, 1601);
    plot(time, IT3);
    title('x[n]');
    % Σχεδιασμός x[n] - x[n-1]
    subplot(2, 1, 2);
    plot(time, FirstDif);
    title('x[n]-x[n-1]');
    
    xn_audio = audioplayer(IT3, 8000); 
    difference_audio = audioplayer(FirstDif, 8000); 
    pause(4);
    play(xn_audio);
    pause(3);
    play(difference_audio);
    audiowrite('x(n).wav', IT3, 8000);
    audiowrite('differenceaudio.wav', FirstDif, 8000);
    FFTFirstDif = fft(FirstDif);
    figure(22);
    plot(freq, abs(FFTFirstDif));
    title('FFT of x[n]-x[n-1]')

% ΆΣΚΗΣΗ 4
    % Ερώτημα α
    cat = imread('image.png');
    figure(23)
    imshow(cat)
    
    % Ερώτημα β
    [x, y] = find(cat,1);
    P = [x, y];
    B = bwtraceboundary(cat,P,'n');
    figure(24)
    plot(B);
    
    % Ερώτημα γ
    x = B(:, 1); 
    y = B(:, 2); 
    
    z = x + 1i*y;
    Z = fft(z);
    figure(25)
    plot(abs(Z));
    
    % Ερώτημα δ
    N=length(Z);
    M1=10;
    M2=50;
    M3=200;
    
    %Για Μ=10
    for n=1:N
        zm(n) = Z(1);
        for k=2:M1
            zm(n)=zm(n)+((1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N));
        end
    end
    
    xm=real(zm);
    ym=imag(zm);
    figure(26)
    subplot(3,2,1);
    plot(xm, ym);
    title('M=10');

    %Για Μ=50
    for n=1:N
        zm(n) = Z(1);
        for k=2:M2
            zm(n)=zm(n)+((1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N));
        end
    end

    xm=real(zm);
    ym=imag(zm);
    subplot(3,2,2);
    plot(xm, ym);
    title('M=50');

    %Για Μ=200
    for n=1:N
        zm(n) = Z(1);
        for k=2:M3
            zm(n)=zm(n)+((1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N));
        end
    end

    xm=real(zm);
    ym=imag(zm);
    subplot(3,2,3);
    plot(xm, ym);
    title('M=200');
    
    % Ερώτημα ε
    N=length(Z);
    M1=10;
    M2=50;
    M3=200;
    count1 = 0; count2 = 0; count3 = 0;

    %Θα δοκιμάσουμε διαφορετικές τιμές για το Μ, προκειμένου να
    %αποφανθούμε ποια από αυτές είναι επαρκής
    %Μ=10
    for n=1:N
        zhatm(n)=Z(1); 
        for k=2:M1/2
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        for k=N-M1/2:N
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        count1 =  count1 +1
    end
    
    yhatm = real(zhatm);
    xhatm = imag(zhatm);
    figure(27)
    subplot(3,2,1)
    plot(xhatm, -yhatm);
    title('M=10');

    %M=50
    for n=1:N
        zhatm(n)=Z(1); 
        for k=2:M2/2
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        for k=N-M2/2:N
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        count2 =  count2 +1
    end
    
    yhatm = real(zhatm);
    xhatm = imag(zhatm);
    subplot(3,2,2)
    plot(xhatm, -yhatm);
    title('M=50');

    %M=200
    for n=1:N
        zhatm(n)=Z(1); 
        for k=2:M3/2
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        for k=N-M3/2:N
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        count3 =  count3 +1
    end
    
    yhatm = real(zhatm);
    xhatm = imag(zhatm);
    subplot(3,2,3)
    plot(xhatm, -yhatm);
    title('M=200');

    %Ερώτημα στ : επαναλαμβάνουμε όλα τα ερωτήματα βήμα βήμα.
    % Ερώτημα στ.α)
    europe = imread('image 2 - refined.png');
    figure(28)
    imshow(europe)
    
    % Ερώτημα στ.β)
    europe_binary = im2bw(europe, 0); %μετατροπή της εικόνας 'image 2 -refined.png' σε binary εικόνα
    [x, y] = find(europe_binary,1);
    figure(29)
    subplot(1,2,1), imshow(europe_binary)
    P = [x, y];
    B = bwtraceboundary(europe_binary,P,'n');
    subplot(1,2,2), plot(B);
    
    % Ερώτημα στ.γ)
    x = B(:, 1); 
    y = B(:, 2); 
    
    z = x + 1i*y;
    Z = fft(z);
    figure(30)
    plot(abs(Z));
    
    % Ερώτημα στ.δ)
    N=length(Z);
    M1=10;
    M2=50;
    M3=200;
    
    %Για Μ=10
    for n=1:N
        zm(n) = Z(1);
        for k=2:M1
            zm(n)=zm(n)+((1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N));
        end
    end
    
    xm=real(zm);
    ym=imag(zm);
    figure(31)
    subplot(3,2,1);
    plot(xm, ym);
    title('M=10');

    %Για Μ=50
    for n=1:N
        zm(n) = Z(1);
        for k=2:M2
            zm(n)=zm(n)+((1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N));
        end
    end

    xm=real(zm);
    ym=imag(zm);
    subplot(3,2,2);
    plot(xm, ym);
    title('M=50');

    %Για Μ=200
    for n=1:N
        zm(n) = Z(1);
        for k=2:M3
            zm(n)=zm(n)+((1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N));
        end
    end

    xm=real(zm);
    ym=imag(zm);
    subplot(3,2,3);
    plot(xm, ym);
    title('M=200');
    
    % Ερώτημα στ.ε)
    N=length(Z);
    M1=10;
    M2=50;
    M3=200;
    count1 = 0; count2 = 0; count3 = 0;

    %Θα δοκιμάσουμε διαφορετικές τιμές για το Μ, προκειμένου να
    %αποφανθούμε ποια από αυτές είναι επαρκής
    %Μ=10
    for n=1:N
        zhatm(n)=Z(1); 
        for k=2:M1/2
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        for k=N-M1/2:N
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        count1 =  count1 +1
    end
    
    yhatm = real(zhatm);
    xhatm = imag(zhatm);
    figure(32)
    subplot(3,2,1)
    plot(xhatm, -yhatm);
    title('M=10');

    %M=50
    for n=1:N
        zhatm(n)=Z(1); 
        for k=2:M2/2
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        for k=N-M2/2:N
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        count2 =  count2 +1
    end
    
    yhatm = real(zhatm);
    xhatm = imag(zhatm);
    subplot(3,2,2)
    plot(xhatm, -yhatm);
    title('M=50');

    %M=200
    for n=1:N
        zhatm(n)=Z(1); 
        for k=2:M3/2
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        for k=N-M3/2:N
            zhatm(n)=zhatm(n)+(1/N)*Z(k).*exp(j*2*pi*(k-1)*n/N);
        end
        count3 =  count3 +1
    end
    
    yhatm = real(zhatm);
    xhatm = imag(zhatm);
    subplot(3,2,3)
    plot(xhatm, -yhatm);
    title('M=200');

    figure(33) 
    europe_perimeter = bwperim(europe_binary)
    imshow(europe_perimeter)