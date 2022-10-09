clear,clc; close all;

% Settings
PER_UNIT_MASS = false;

% Include subfolders
addpath('data\');
addpath('sections\');
addpath('utils\');
addpath('sulpis2021\');

% Color schemes
DEFAULT_COLORS = [0, .447, .741; .85, .325, .098; .929, .694, .125; .494, .184, .556; .466, .674, .188; .301, .745, .933; .635, .078, .184];
FLUX_COLORS = (1./256) * [
    49,130,189; % E_bdy
    222,45,38;  % E_meso
    252,146,114;% E_iso
    49,163,84;  % M_meso
    161,217,155;% M_iso
    82,82,82;   % total WTF
    ];
BIOME_COLORS = load('biome_colors.mat').colors;
BLUES = [[189,215,231];[107,174,214];[49,130,189];[8,81,156]] / 255;
GREENS = [[186,228,179];[116,196,118];[49,163,84];[0,109,44]] / 255;
REDS = [[252,174,145];[251,106,74];[222,45,38];[165,15,21]] / 255;
PURPLES = [[203,201,226];[158,154,200];[117,107,177];[84,39,143]] / 255;
ORANGES = [[253,190,133];[253,141,60];[230,85,13];[166,54,3]] / 255;

% Other constants
sec_in_year = 365.25 * 24 * 3600;
ALK = 1; SAL = 2; NIT = 3; ALKP = 4; ALKC_lin = 5; ALKC_reg = 6; ALKC_avg = 7; % enum for easy indexing
save('data\constants.mat');


%% Load data
[sulpis, wmt_data, wmt_alkP_kd] = gather_data();
data = process_data(sulpis, wmt_data, wmt_alkP_kd);


%% Show figures
tic;

section21_seawater_transformation(data);
section22_biome_map();
section24_density_depth_conversion();
section31_dissolution_profiles_alkstar(data);
section32_dissolution_profiles_wmt(data);
section33_error_compensation(data);
section34_salinity_normalization(data);
section35_diffusivity_uncertainty(data);

toc;



%% FUNCTIONS
function [sulpis, wmt_data, wmt_alkP_kd] = gather_data()
    load('constants.mat');
    % Sulpis 2021
    if PER_UNIT_MASS
        sulpis_fname = "processed_sulpis_perUnitMass.mat";
    else
        sulpis_fname = "processed_sulpis.mat";
    end
    if exist(sulpis_fname, "file")
        load(sulpis_fname, 'sulpis');
    else
        sulpis_raw = load('D:\NIOZ\Sulpis_reproduction\pelagic_dissolution-main\Sulpis2021_data_gamma_withNegatives_extendedBiomes.mat');
        sulpis = process_sulpis(sulpis_raw);
        save("sulpis2021\"+sulpis_fname, 'sulpis');
    end
    
    % WMT
    wmt_data = struct('ndensity', [], 'Grho', [], 'Gtracer', [], 'dia_M', []);
    suffix = ["", "_pureSWR"];
    folder = "data\WMT_results\";
    for i = 1:2
        wmt_data(ALK, i) = load(folder + "WMT_results_Alk_varK"+suffix(i)+".mat");
        wmt_data(SAL, i) = load(folder + "WMT_results_SA_varK"+suffix(i)+".mat");
        wmt_data(NIT, i) = load(folder + "WMT_results_Nit_varK"+suffix(i)+".mat");
        wmt_data(ALKC_lin, i) = load(folder + "WMT_results_AlkC_linear_varK"+suffix(i)+".mat");
        wmt_data(ALKC_reg, i) = load(folder + "WMT_results_AlkC_biomeFits_varK"+suffix(i)+".mat");
        wmt_data(ALKC_avg, i) = load(folder + "WMT_results_AlkC_average_varK"+suffix(i)+".mat");
    end

    % WMT of AlkP for varying eddy diffusivities K and D
    wmt_alkP_kd = struct('ndensity', [], 'Grho', [], 'Gtracer', [], 'dia_M', []);
    for k = 1:3
        for d = 1:3
            kval = 0.5 * k;
            dval = 0.5 * d;
            wmt_alkP_kd(k,d) = load(folder + "WMT_results_AlkP_"+sprintf("%0.2f", kval)+"K"+sprintf("%0.2f", dval)+"D.mat");

            % WMT error compensation
            wmt_sal = load(folder + "WMT_results_SA_"+sprintf("%0.2f", kval)+"K"+sprintf("%0.2f", dval)+"D.mat");
            wmt_alkP_kd(k,d).Gtracer = wmt_alkP_kd(k,d).Gtracer - 66.1 * wmt_sal.Gtracer;
            wmt_alkP_kd(k,d).dia_M = wmt_alkP_kd(k,d).dia_M - 66.1 * wmt_sal.dia_M;
        end
    end
end

function [data] = process_data(sulpis, wmt_data, wmt_alkP_kd)
    load('constants.mat');
    [wmt_data] = process_wmt(wmt_data);
    [accum_data] = compute_accumulation(wmt_data);
    [wmt_alkP_kd] = process_wmt(wmt_alkP_kd);
    [accum_alkP_kd] = compute_accumulation(wmt_alkP_kd);
    
    for i = 1:size(wmt_data, 2)
        wmt_data(ALKP,i).gamma = wmt_data(ALK,i).gamma;
        wmt_data(ALKP,i).E = wmt_data(ALK, i).E + 1.26 * wmt_data(NIT, i).E;
        wmt_data(ALKP,i).M = wmt_data(ALK, i).M + 1.26 * wmt_data(NIT, i).M;
        wmt_data(ALKP,i).tot = wmt_data(ALK, i).tot + 1.26 * wmt_data(NIT, i).tot;
        accum_data(ALKP,i).gamma = accum_data.gamma;
        accum_data(ALKP,i).E = accum_data(ALK, i).E + 1.26 * accum_data(NIT, i).E;
        accum_data(ALKP,i).M = accum_data(ALK, i).M + 1.26 * accum_data(NIT, i).M;
        accum_data(ALKP,i).tot = accum_data(ALK, i).tot + 1.26 * accum_data(NIT, i).tot;
    end
    data = struct('wmt_data', wmt_data, 'wmt_alkP_kd', wmt_alkP_kd, 'accum_data', accum_data, 'accum_alkP_kd', accum_alkP_kd, 'sulpis', sulpis);
end


function [sulpis] = process_sulpis(raw_struct)
    load('constants.mat');
    disp("Processing raw Sulpis data...");
    dgamma = 0.05;
    gamma_edges = 15:dgamma:30;
    gamma_array = gamma_edges(1:end-1) + dgamma / 2;
    sulpis = struct('gamma', gamma_array, 'd_rate', nan(size(gamma_array)), 'd_rate_std', nan(size(gamma_array)), 'd_rate_biomes', nan([10, length(gamma_array)]), 'd_rate_biomes_std', nan([10, length(gamma_array)]));
    
    for i = 1:length(gamma_array)
        I = raw_struct.Gamma_mean > gamma_edges(i) & raw_struct.Gamma_mean <= gamma_edges(i+1);
        sulpis.d_rate(i) = mean(raw_struct.Alkstar_rate(I), 'omitnan');
        sulpis.d_rate_std(i) = std(raw_struct.Alkstar_rate(I), 'omitnan');
    end
    
    for j = 1:10
        gmj = squeeze(raw_struct.Gamma_mean(:,j,:));
        arj = squeeze(raw_struct.Alkstar_rate(:,j,:));
        for i = 1:length(gamma_array)
            I = gmj > gamma_edges(i) & gmj <= gamma_edges(i+1);
            sulpis.d_rate_biomes(j,i) = mean(arj(I), 'omitnan');
            sulpis.d_rate_biomes_std(j,i) = std(arj(I), 'omitnan');
        end
    end
    
    % Smoothen curves
    wd_sz = 3; % number of points on either side
    disp("Sulpis: "+sprintf("%d", wd_sz*2+1)+" point average");
    sulpis.d_rate = moving_average(sulpis.d_rate, wd_sz);
    sulpis.d_rate_std = moving_average(sulpis.d_rate_std, wd_sz);
    sulpis.d_rate_biomes = moving_average(sulpis.d_rate_biomes, wd_sz);
    sulpis.d_rate_biomes_std = moving_average(sulpis.d_rate_biomes_std, wd_sz);
    
    idx = find(~isnan(sulpis.d_rate), 1);
    sulpis.d_rate(idx) = nan; % remove first point, it only originates from Indian Ocean data and is therefore an outlier

    % Remove NaN values
    I_not_nan = ~isnan(sulpis.d_rate + sulpis.d_rate_std);
    sulpis.gamma = sulpis.gamma(I_not_nan);
    sulpis.d_rate_std = sulpis.d_rate_std(I_not_nan);
    sulpis.d_rate = sulpis.d_rate(I_not_nan);
    sulpis.d_rate_biomes = sulpis.d_rate_biomes(:,I_not_nan);
    sulpis.d_rate_biomes_std = sulpis.d_rate_biomes_std(:,I_not_nan);

    sulpis.mass = mass_for_densities(sulpis.gamma);
    sulpis.biomes_mass = mass_for_densities_ext_regions(sulpis.gamma);
    if ~PER_UNIT_MASS
        % umol/kg/yr -> Tmol/(kg/m^-3)/yr
        sulpis.d_rate = sulpis.d_rate .* sulpis.mass * 1e-18 / dgamma;
        sulpis.d_rate_std = sulpis.d_rate_std .* sulpis.mass * 1e-18 / dgamma;
        sulpis.d_rate_biomes = sulpis.d_rate_biomes .* sulpis.biomes_mass * 1e-18 / dgamma;
        sulpis.d_rate_biomes_std = sulpis.d_rate_biomes_std .* sulpis.biomes_mass * 1e-18 / dgamma;
    end
    disp('done!');
end


function [wmt] = process_wmt(wmt)
    load('constants.mat');
    for i = 1:size(wmt, 1)
        for j = 1:size(wmt, 2)
            % Flip sign (we define positive value as WMT towards higher densities)
            wmt(i, j).gamma = wmt(i, j).ndensity;
            wmt(i, j).E = -wmt(i, j).Gtracer;
            wmt(i, j).M = -wmt(i, j).dia_M;
            wmt(i, j).G = -wmt(i, j).Grho;

            wmt(i, j).E(1, wmt(i,j).gamma < 16) = 0;  % noise at impossible density
       
            % umol/s -> Tmol/yr
            wmt(i, j).E = sec_in_year * 1e-18 * wmt(i, j).E;
            wmt(i, j).M = sec_in_year * 1e-18 * wmt(i, j).M;
            
            % Total
            wmt(i, j).tot = sum(wmt(i, j).E, 1) + sum(wmt(i, j).M, 1);
        end
    end
end

function [accum] = compute_accumulation(wmt)
    load('constants.mat');
    accum = struct([]);
%     mass = mass_for_densities(wmt(1,1).gamma(1:end-1) + diff(wmt(1,1).gamma));
%     save('mass_isopycnals.mat', 'mass');
    load('mass_isopycnals.mat', 'mass');
    for i = 1:size(wmt, 1)
        for j = 1:size(wmt, 2)
            % Compute accumulation
            accum(i, j).gamma = wmt(i, j).gamma(1:end-1) + diff(wmt(i, j).gamma);
            accum(i, j).E = -diff(wmt(i, j).E, 1, 2); % E(x - x/2) - E(x + x/2)
            accum(i, j).M = -diff(wmt(i, j).M, 1, 2); % M(x - x/2) - M(x + x/2)
            accum(i, j).G = -diff(wmt(i, j).G, 1, 2); % G(x - x/2) - G(x + x/2)
            
            if PER_UNIT_MASS && ~isempty(accum(i, j).E)
                % Tmol/yr -> umol/kg/yr
                accum(i, j).E = 1e18 * accum(i, j).E ./ mass;
                accum(i, j).M = 1e18 * accum(i, j).M ./ mass;

                accum(i,j).E(:, accum(i,j).gamma > 28.4) = 0;
                accum(i,j).M(:, accum(i,j).gamma > 28.4) = 0;
            else
                % Tmol/yr -> Tmol/(kg m^-3)/yr
                accum(i, j).E = accum(i, j).E ./ diff(wmt(i,j).gamma);
                accum(i, j).M = accum(i, j).M ./ diff(wmt(i,j).gamma);
            end

            % Compute totals
            accum(i, j).tot = sum(accum(i, j).E, 1) + sum(accum(i, j).M, 1);
        end
    end
end



