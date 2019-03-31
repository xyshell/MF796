function myplot2(pdf, k, T)
% plot risk neutral density vs strike
% pdf: risk neutral density
% k: strike
pdfK = k(3:end);
for i = 1:size(pdf, 2)
    subplot(1, 2, i)
    plot(pdfK, pdf(:, i), 'LineWidth', 2)
    xlabel('Strike (K)')
    ylabel('Value')
    title(['T = ', num2str(T(i))])
    grid
end
end

