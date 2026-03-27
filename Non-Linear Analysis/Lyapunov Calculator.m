t = RM1.Time;
%%
cut_time = 0.5;
idx = t > cut_time;

t_ss = t(idx);
Vc_ss = Vc.Data(idx);
Vc1_ss = Vc1.Data(idx);
Vc2_ss = Vc2.Data(idx);
Vc3_ss = Vc3.Data(idx);
Rm1_ss = RM1.Data(idx);
Rm2_ss = RM2.Data(idx);
u_ss = Input.Data(idx);

[~,poincare_idx] = findpeaks(u_ss);

figure('Name', 'Dynamical System Phase Portraits', 'NumberTitle', 'off');

% 1. The 2D Voltage Phase Portrait (Classic view)
subplot(2,2,1);
plot(Vc_ss, Vc1_ss, 'b', 'LineWidth', 0.5);
xlabel('Voltage C (V)');
ylabel('Voltage C2 (V)');
title('2D State Space: V_{C} vs V_{C1}');
grid on;

% 2. The 3D Voltage Attractor (Best for spotting chaos)
subplot(2,2,2);
plot3(Vc_ss, Vc1_ss, Vc3_ss, 'r', 'LineWidth', 0.5);
xlabel('Voltage C (V)');
ylabel('Voltage C1 (V)');
zlabel('Voltage C3 (V)');
title('3D Attractor');
grid on;
view(45, 30); % Adjust viewing angle

% 3. The Memristor State Space (Voltage vs Resistance)
subplot(2,2,3);
plot(Vc1_ss - Vc_ss, Rm1_ss, 'k', 'LineWidth', 0.5);
xlabel('Voltage across Memristor 1 (V)');
ylabel('Memristance R_{M1} (\Omega)');
title('Memristor Pinched Hysteresis');
grid on;

% 3. The Memristor 2 State Space (Voltage vs Resistance)
subplot(2,2,4);
plot(Vc_ss - Vc3_ss , Rm2_ss, 'k', 'LineWidth', 0.5);
xlabel('Voltage across Memristor 2 (V)');
ylabel('Memristance R_{M2} (\Omega)');
title('Memristor Pinched Hysteresis');
grid on;

Vc_poincare  = Vc_ss(poincare_idx);
Vc1_poincare = Vc1_ss(poincare_idx);

figure('Name', 'Poincaré Section', 'Color', 'w');
plot(Vc_poincare, Vc1_poincare, '.', 'MarkerSize', 20, 'Color', [0 0.4470 0.7410]); 
xlabel('Voltage C (V)');
ylabel('Voltage C1 (V)');
title('Poincaré Section (Sampled at Input Drive Period)');
grid on;


%%
model_name = 'Adaptive_Oscillator.slx'; % Replace with your actual .slx name
% --- Configuration ---

A_sweep = linspace(0, 1.0, 100); % Sweep amplitude from 0.1V to 1.5V
% ---------------------

figure('Name', 'Bifurcation Diagram', 'Color', 'w');
hold on;
xlabel('Sine Wave Amplitude (V)');
ylabel('Peaks of V_C (V)');
title('Bifurcation Diagram: Route to Chaos');
grid on;

disp('Starting parameter sweep... this may take a few minutes.');

for i = 1:length(A_sweep)
    % 1. Set the variable in the base workspace so Simulink sees it
    A_in = A_sweep(i);
    disp(A_in);
    % 2. Run the simulation silently
    out = sim(model_name, 'SimulationMode', 'accelerator', 'ReturnWorkspaceOutputs', 'on');
    
    % 3. Extract data robustly (adjust 'Vc' to match your workspace variable name)
    % Assuming the data is in out.Vc as a timeseries object:
    time_vec = out.tout;
    vc_raw = out.logsout{4}.Values.Data;
    
    % 4. Strip the transient
    idx = time_vec > 0; 
    vc_ss = vc_raw(idx);
    
    % 5. Find the peaks of the steady-state waveform
    peaks = findpeaks(vc_ss);
    
    % 6. Plot the peaks for this specific amplitude
    % We plot them as tiny black dots. If there are 2 peaks, we get 2 dots. 
    % If it's chaotic, we get a vertical spread of dots.
    plot(A_sweep(i) * ones(size(peaks)), peaks, 'k.', 'MarkerSize', 25);
    
    % Optional: Force MATLAB to draw the plot as it calculates
    drawnow; 
end

disp('Sweep complete.');
hold off;
%%
% --- 1. Setup Data ---
% Make sure vc_ss is a column vector
data = vc_ss(:) - mean(vc_ss); 
dt = mean(diff(time_vec));
N_total = length(data);

% --- 2. Find Time Delay (tau) via Autocorrelation ---
disp('Calculating time delay...');
max_lag = 500;
acf = zeros(max_lag, 1);
var_data = sum(data.^2);

for lag = 1:max_lag
    % Manual autocorrelation (no toolboxes needed)
    acf(lag) = sum(data(1:end-lag) .* data(lag+1:end)) / var_data;
end

% Find where autocorrelation drops below 1/e (~0.36)
tau = find(acf < exp(-1), 1);
if isempty(tau)
    tau = 15; % Fallback heuristic if it oscillates slowly
end
fprintf('Optimal time delay (tau) = %d samples\n', tau);

% --- 3. Reconstruct Phase Space ---
m = 5; % Embedding dimension (matches your circuit's complexity)
N = N_total - (m-1)*tau;
Y = zeros(N, m);

for i = 1:m
    Y(:, i) = data((1:N) + (i-1)*tau);
end

% --- 4. Track Divergence (Rosenstein Algorithm) ---
disp('Tracking trajectory divergence... (This may take a minute)');
track_steps = 500; % How many steps into the future to track
d_m = zeros(track_steps, 1);
valid_points = 0;

theiler_window = 100; % Ignore points from the same immediate orbit
num_eval_points = 500; % Evaluate 500 reference points to save computation time
step_size = floor((N - track_steps) / num_eval_points);

for i = 1:step_size:(N - track_steps)
    ref_pt = Y(i, :);
    
    % Calculate squared distances to all other points
    diffs = Y - ref_pt;
    sq_dists = sum(diffs.^2, 2);
    
    % Apply Theiler window (set local temporal points to infinity)
    exclude_start = max(1, i - theiler_window);
    exclude_end   = min(N, i + theiler_window);
    sq_dists(exclude_start:exclude_end) = inf;
    
    % Find the nearest neighbor in the phase space
    [~, nn_idx] = min(sq_dists);
    
    % Track how far apart they drift over 'track_steps'
    if nn_idx + track_steps <= N
        for k = 1:track_steps
            dist_k = norm(Y(i+k, :) - Y(nn_idx+k, :));
            if dist_k > 0
                d_m(k) = d_m(k) + log(dist_k);
            end
        end
        valid_points = valid_points + 1;
    end
end

% Average the logarithmic divergence
d_m = d_m / valid_points;
time_axis = (1:track_steps) * dt;

% --- 5. Extract Slope (LLE) and Plot ---
% We calculate the slope of the initial linear rise. 
% You may need to adjust 'fit_end' based on the plot.
fit_start = 5;
fit_end = floor(track_steps * 0.3); % Fit the first 30% of the curve

p = polyfit(time_axis(fit_start:fit_end)', d_m(fit_start:fit_end), 1);
LLE = p(1);

figure('Name', 'Lyapunov Exponent', 'Color', 'w');
plot(time_axis, d_m, 'k', 'LineWidth', 1.5); hold on;
plot(time_axis(fit_start:fit_end), polyval(p, time_axis(fit_start:fit_end)), 'r--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Average Log Divergence <ln(d)>');
title(sprintf('Largest Lyapunov Exponent (LLE) = %.4f', LLE));
grid on; hold off;

fprintf('====================================\n');
fprintf('Largest Lyapunov Exponent: %.4f\n', LLE);
fprintf('====================================\n');