# This is the readme with function descriptions

## File DATA_generate_simulation.R

Function **generate_simulation_data**
- generates simulation data from equation 5 of write up

## File FUNC_paramater_estimates.R

Function **clean_y**
- removes NA values from data
- necessary for removing NAs from data with days where flight did not occur

Function **get_matern_values**
- function to generate Matern(5/2) values

Function **get_matern**
- calls *get_matern_values* to get the actual Matern kernel

Function **get_K_i**
- returns the K matrix, sigma_2 * M from *get_matern*

Function **get_V_i**
- returns the V matrix
- calls *get_K_i* and adds sigma_2 * identity matrix


Function **get_V_i**
- returns the V matrix
- calls *get_K_i* and adds sigma_2 * identity amtrix

Function **get_h_j**
- calculates h vector from equation 6
- TODO this may have some issues 

Function **get_g_i**
- calculates g values from equation 6
- TODO this may have some issues 

Function **get_g**
- calculates actual g vector from equation 6
- calls *get_g_i* to get values
- TODO this may have some issues

Function **get_sigma_mu_post**
- function to calculate sigma_mu_post from full conditional of mu_i

Function **get_alpha_mu_post**
- TODO may have sigma_2^(-2) when we really just need sigma_2(^-1)
- this may be wrong, going to have someone look at it

Function **get_mu_i**
- function to draw mu values
- calls *get_alpha_mu_post* and *get_sigma_mu_post* to get the mean and variance

Function **vector_differences**
- this is the function to calculate Yi = mu_i - gI
- used in calculating the sigma_2 updates

Function **get_sigma_squared**
- function to do IG draw for sigma_2
- calls *vector_differences* to get one of the terms needed

Function **lk_acceptance**
- function to get alpha(lk, lk')

Function **get_lk**
- function to do MH update of lk
- calls *lk_acceptance* to get values for udpates

Function **get_H_matrix**
- function to create H matrix for the xi draws
- calls *get_h_j* for individual values

Function **psi_xi**
- function to calculate negative log likelihood for xi values

Function **get_xi**
- function to update the xi values
- calls *psi_xi* for acceptance values
- calls *samp.WC* (Pallavi's code) for initial xi values

Function **psi_alpha**
- function for calculating likelihood of alpha

Function **get_alpha**
- function for calculating alpha values
- calls **psi_alpha** for updates

Function **lb_acceptance**
- function to get alpha(lb, lb')

Function **get_lb**
- function to do MH update of lb
- calls *lb_acceptance* to get values for udpates




