function simulate_amplifier()

    % Simulation constants
    SAMPLE_RATE = 44100;
    SIGNAL_LENGTH = 44100;
    SIGNAL_FREQ_1 = 20;
    SIGNAL_FREQ_2 = 2000;
    SIGNAL_FREQ_3 = 4;
    SIGNAL_AMPLITUDE_1 = 1;
    SIGNAL_AMPLITUDE_2 = 0.3;
    SIGNAL_AMPLITUDE_3 = 1;

    % Input signal as a series of pulses
    input_signal = ones(1, SIGNAL_LENGTH);
    input_signal(10:20) = fliplr(0 : 0.1 : 1);
    input_signal = (SIGNAL_AMPLITUDE_1 * sin(2 * pi * SIGNAL_FREQ_1 * [0 : SIGNAL_LENGTH-1] / SAMPLE_RATE) + ...
                   SIGNAL_AMPLITUDE_2 * sin(2 * pi * SIGNAL_FREQ_2 * [0 : SIGNAL_LENGTH-1] / SAMPLE_RATE)) .* ...
                   (SIGNAL_AMPLITUDE_3 * sin(2 * pi * SIGNAL_FREQ_3 * [0 : SIGNAL_LENGTH-1] / SAMPLE_RATE));

    % Spectrum of input signal
    fft_in = (fft(input_signal));
    % Process input signal with the filter
    fft_h = H_triode_classA([(0 : SIGNAL_LENGTH/2)]);
    fft_h = [fft_h conj(fliplr(fft_h(2:length(fft_h) -1)))];
    % Obtain output signal
    fft_out = fft_in .* fft_h;
    output_signal = real(ifft(fft_out));
    

    
%     subplot(2,1,1)
%     plot(phase(fft_h))
%     title('Amp');
%     xlabel('Phase')
%     grid
%     subplot(2,1,2)
%     plot(20*log10(abs(fft_h)))
%     xlabel('Magnitude')
%     grid
%     return

    
%     subplot(2,1,1)
%     plot(input_signal(1:50));
%     grid
%     subplot(2,1,2)
%     plot(output_signal(1:50));
%     grid
% return

    plot_handler = figure();
    set(gcf, 'Position', get(0,'Screensize'));
    frequency_axis = [0 : SAMPLE_RATE / SIGNAL_LENGTH : SAMPLE_RATE-1];

    subplot(2,3,1)
    plot(input_signal)
    ylim([-2 2])
    xlim([0 2000])
    grid
    title('Input signal')

    subplot(2,3,2)
    semilogx(frequency_axis, 20*log10(abs(fft_h)))%, ':', 'LineWidth', 2)
    ylim([28 34])
%     xlim([1 SAMPLE_RATE])
    grid
    title('Filter response - Module')
    subplot(2,3,5)
    semilogx(frequency_axis, unwrap(phase(fft_h)) * 180 / pi)%, ':', 'LineWidth', 2)
%     ylim([32 45])
%     xlim([1 SAMPLE_RATE])
    grid
    title('Filter response - Phase')

    
    
    subplot(2,3,3)
    plot(output_signal)
    ylim([-100 100])
    xlim([0 2000])
    grid
    title('Output signal')

    subplot(2,3,4)
    semilogx(frequency_axis, 20*log10(abs(fft_in)))%, 'x-')
    ylim([-300 100])
    grid
    title('Spectrum input')

    subplot(2,3,6)
    semilogx(frequency_axis, 20*log10(abs(fft_out)))%, 'x-')
    ylim([-300 100])
    grid
    title('Spectrum output')

end
