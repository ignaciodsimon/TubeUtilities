function filter_method_proof()
    % This script was made to prove that implementing the response of the
    % complete circuit for any random input signal can be done by
    % considering the input-output transfer function as the coefficients of
    % the FFT of the filter that needs to be applied to the input signal.
    %
    % To test it, the input was a signal with a DC component a series of
    % impulses. The implemented transfer function was a simple high pass RC
    % filter. The output does not have DC and shows the transitory effects
    % very well.
    %
    % Joe.

    SAMPLE_RATE = 44100;
    SIGNAL_LENGTH = 1000;
    SIGNAL_FREQ = 1000;

    % Input signal as a sine
    input_signal = sin(2 * pi * SIGNAL_FREQ * [0 : SIGNAL_LENGTH-1] / SAMPLE_RATE);

    % Input signal as a series of pulses
    input_signal = ones(1, SIGNAL_LENGTH);
    input_signal(1) = 0;
    input_signal(10) = 0;
    input_signal(12) = 0;
    input_signal(17) = 0;

    % Process input signal with the filter
    fft_in = fft(input_signal);
    fft_h = H([0 : SIGNAL_LENGTH-1]);
    fft_out = fft_in .* fft_h;

    % Obtain output signal
    output_signal = real(ifft(fft_out));

    
    subplot(2,3,1)
    plot(input_signal)
    ylim([-2 2])
    grid
    title('Input signal')

    subplot(2,3,2)
    semilogx(20*log10(abs(fft_h)))
    grid
    title('Filter response')

    subplot(2,3,3)
    plot(output_signal)
    ylim([-2 2])
    grid
    title('Output signal')

    subplot(2,3,4)
    semilogx(20*log10(abs(fft_in)))
    grid
    title('Spectrum input')

    subplot(2,3,6)
    semilogx(20*log10(abs(fft_out)))
    grid
    title('Spectrum output')

end

function weights = H(frequencies)
    % Weights (FFT) of a simple first-order RC filter.
    C = 10^-6;
    R = 1000;
    
    % Hi-pass:
    weights = 1 ./ (1 + ( 1./(1i * 2 * pi * frequencies * C * R) ));

    % Low-pass:
    %weights = 1 ./ (1 + (1i * 2 * pi * frequencies * C * R) );
end
