from Crypto import Random
from Crypto.Random import random
from Crypto.Util import number
from Crypto.PublicKey import ElGamal
from functools import reduce
import numpy as np 
import secrets
from timeit import default_timer as timer

secretsGenerator = secrets.SystemRandom()

class Verifier:
    
    def __init__(self) -> None:
        self.verifierTime = []
    def generateChallenge(self, params):
        # Get challenge x that is an element of Z_q*
        return secretsGenerator.randint(1, params[0] - 1)
    
    def zero_arg_verify(self, params, com_a0, com_bm, vec_com_D, vec_a, vec_b, r, s, t, vec_x_m, vec_x_m_j, vec_x_m_k, vec_comA, vec_comB, y ):
        start = timer()
        print(f'Is Com_A0 an element of G? : {isComGroupElement(params, com_a0)}')
        print(f'Is Com_Bm an element of G? : {isComGroupElement(params, com_bm)}')
        print(f'Does Com_D contain 2m+1 ({2*m+1}) elements? : {len(vec_com_D)}')
        print(f'Is each component of com_D a valid element of G? : {isGenComGroupElement(params, vec_com_D)}')
        print(f'Is Com_D_(m+1) == com(0;0) ?: {vec_com_D[m+1] == 1}')
        print(f'Does vector a and vector b contain n elements of Z_q?: {len(vec_a) == n and len(vec_b) == n}')

        # Now equations
        vec_comA.insert(0, com_a0)
        comA_power_x = bilinearMapCom(vec_comA, vec_x_m, p)
        com_a_r = pc.commit(vec_a, params, r)
        print(f'Is the product of C_A_i raised to the power x^i == com(a;r) ? : {comA_power_x == com_a_r} ')


        vec_comB.insert(len(vec_comB), com_bm)
        comB_power_x = bilinearMapCom(vec_comB, vec_x_m_j, p)
        com_b_s = pc.commit(vec_b, params, s)
        print(f'Is the product of C_B_j raised to the power x^(m-j) == com(b;s) ? : {comB_power_x == com_b_s} ')

        comD_power_x = bilinearMapCom(vec_com_D, vec_x_m_k, p)
        com_a_cross_b_t = pc.commit(zeroArgBilinearMap(vec_a, vec_b, y), params, t)
        print(f'Is the product of C_D_k raised to the power x^k == com(a * b;t) ? : {comD_power_x == com_a_cross_b_t}')
        end = timer()
        self.verifierTime.append(end-start)

        # To raise an exception if one of the tests fail
        if (len(vec_com_D) == 2*m+1 and len(vec_a) == n and len(vec_b) == n and comA_power_x == com_a_r and comB_power_x == com_b_s and comD_power_x == com_a_cross_b_t) == False:
            raise Exception("Error in zero argument verification")

    def hadamard_arg_verify(self, params, vec_comA, com_b, vec_com_B ):
        start = timer()
        # Check cB2, . . . , cBm−1 ∈ G, cB1 = cA1, cBm = cb
        print(f'Are C_B2,...C_Bm-1 valid elements of G? : {isGenComGroupElement(params, vec_com_B[1:-1])}')
        print(f'Is C_B1 == C_A1 ? : {vec_com_B[0] == vec_comA[0]}')
        print(f'Is C_Bm == C_b ? : {vec_com_B[-1] == com_b}')
        end = timer()
        self.verifierTime.append(end-start)

        # To raise an exception if one of the tests fail
        if (vec_com_B[0] == vec_comA[0] and vec_com_B[-1] == com_b) == False:
            raise Exception("Error in Hadamard argument verification")

    def single_value_product_arg_verify(self, pc, params, comA, b, com_d, com_delta, com_triangle, x_5_3, vec_a_bar, vec_b_bar, r_bar, s_bar):
        start = timer()
        q, _, _, p = params 
        print(f'Is com_d and element of G?: {isComGroupElement(params, com_d)}')
    
        print(f'Is com_delta and element of G?: {isComGroupElement(params, com_delta)}')

        print(f'Is com_triangle and element of G?: {isComGroupElement(params, com_triangle)}')

        com_a_bar = pc.commit(vec_a_bar, params, r_bar)

        vec_x_b_bar = [ (x_5_3*b_bar_i)%q for b_bar_i in vec_b_bar[1:]]
        vec_b_bar_a_bar = [ (b_i * a_j) % q for b_i, a_j in zip(vec_b_bar[:-1], vec_a_bar[1:])]

        vec_x_b_bar_minus_b_bar_a_bar = [ (x - y)%q for x, y in zip(vec_x_b_bar, vec_b_bar_a_bar)] 
        com_last = pc.commit( vec_x_b_bar_minus_b_bar_a_bar, params, s_bar )

        com_a_power_x = pow(comA, x_5_3, p)

        com_a_power_x_times_com_d = (com_a_power_x * com_d)%p
        print(f'Is com_a_power_x * com_d == com_ck(vec_a_bar; r_bar) ?: {com_a_power_x_times_com_d == com_a_bar} ')
    
        com_triangle_power_x = pow(com_triangle, x_5_3, p)
        com_triangle_power_x_times_com_delta = (com_triangle_power_x * com_delta)%p
        print(f'Is com_triangle_power_x * com_delta == com_ck(vec_x_b_bar - vec_b_bar_a_bar; s_bar) ?: {com_triangle_power_x_times_com_delta == com_last} ')

        print(f'Is b_bar_1 == a_bar_1 ? : {vec_a_bar[0] == vec_b_bar[0]}')
        print(f'Is b_bar_n == x*b ? : { (x_5_3*b)%q == vec_b_bar[-1]}')
        end = timer()
        self.verifierTime.append(end-start)
        
        # To raise an exception if one of the tests fail
        if (com_a_power_x_times_com_d == com_a_bar and com_triangle_power_x_times_com_delta == com_last and vec_a_bar[0] == vec_b_bar[0] and (x_5_3*b)%q == vec_b_bar[-1] ) == False:
            raise Exception("Error in Single value product argument verification")

    def product_arg_verify(self, params, com_b): 
        start = timer()
        print(f'Is com_b a valid group element?: {isComGroupElement(params, com_b)}')
        end = timer()
        self.verifierTime.append(end-start)

    def multi_expo_optim_arg_verify(self, pc, params, m, n, com_A0, com_B_k, vec_E, vec_a, r, b, s, tau, x_multi, singleCiphertext, com_A, shuffled_ciphertexts ):
        start = timer()
        q, _, _, p = params 

        vec_x_prime = [pow(x_multi, i, q) for i in range(1, m+1)]
        temp_vec_x = [pow(x_multi, i, q) for i in range(2*m)]
        mat_c = shuffled_ciphertexts

        print(f'Is commitment to A_0 an element of G? : {isComGroupElement(params, com_A0)} ')
        print(  
        f'Are all commitments B_k an element of group G?: {isGenComGroupElement(params, com_B_k)}')
        print(f'Are all E_k`s elements of H: G x G ?: {groupElement(params, vec_E )}')
        print(f'Does vec_a have n ({n}) components in Z_q?: {len(vec_a) == n}')
        print(f'Are r, b, s and tau elements of Z_q ?: {0 <= r <= q-1 and 0 <= b <= q-1 and 0 <= s <= q-1 and 0<=tau <= q-1 }')
        print(f'Is c_B_m == com_ck(0;0)? : {com_B_k[m] == 1}')


        print(f'Is E_m == C? : {vec_E[m] == singleCiphertext}')

        # Now for equation verifications
        com_a_pow_x = bilinearMapCom(com_A,vec_x_prime, p)
        com_a0_times_a_pow_x = (com_A0 * com_a_pow_x)%p
        com_a_r = pc.commit(vec_a, params, r)
        print(f'Is com_a0 * com_a^x == com_a_r?: {com_a0_times_a_pow_x == com_a_r}')

        com_b_s = pc.commit(b, params, s)
        com_B_k_pow_x = bilinearMapCom(com_B_k,temp_vec_x, p) 
        print(f'Is com(b;s) == com_B_k^(x^k) ?: {com_B_k_pow_x== com_b_s}')

        e_k_power_x = bilinearMap(vec_E, temp_vec_x, p)
        enc_k_temp = elG._encrypt(pow(elG.publickey().g, b, p), tau)
        prod = [1,1]
        temp_vec_x_2 = [pow(x_multi, m-i-1,q) for i in range(m)]

        # Multiply each component of the power of x with each entry of vector a 
        result = [ (a * b)%q for a in temp_vec_x_2 for b in vec_a]
        ''' Result should contain a vector of length n*m. We need to split this 
            into components of n again to get m vectors in total
            Each of which is of size n'''
        pow_rhs = split_list(result, n)

        for i in range(m):
             tmp = bilinearMap(mat_c[i], pow_rhs[i], p)
             prod = [(x*y)%p for x,y in zip(prod, tmp)]
        # prod should contain C_i^(x^(m−i)*a).
        rhs_eq = [(i*j)%p for i,j in zip(enc_k_temp, prod)]
        print(f'Does the product of E_k`s raised to the power x^k == the right hand side ?: {e_k_power_x == rhs_eq}')
        end = timer()
        self.verifierTime.append(end-start)
        if (len(vec_a) == n and 0 <= r <= q and 0 <= b <= q and 0 <= s <= q and 0<=tau <= q and com_B_k[m] == 1 and vec_E[m] == singleCiphertext and com_a0_times_a_pow_x == com_a_r and com_B_k_pow_x== com_b_s and e_k_power_x == rhs_eq) == False:
            raise Exception("Error in Multi exponentiation modified verification")

    def multi_expo_arg_verify(self, pc, params, m, n, com_A0, com_B_k, vec_E, vec_a, r, b, s, tau, x_multi, singleCiphertext, com_A, shuffled_ciphertexts ):
        start = timer()
        q, _, _, p = params 

        vec_x_prime = [pow(x_multi, i, q) for i in range(1, m+1)]
        temp_vec_x = [pow(x_multi, i, q) for i in range(2*m)]
        mat_c = split_list(shuffled_ciphertexts, n )

        print(f'Is commitment to A_0 an element of G? : {isComGroupElement(params, com_A0)} ')
        print(  
        f'Are all commitments B_k an element of group G?: {isGenComGroupElement(params, com_B_k)}')
        print(f'Are all E_k`s elements of H: G x G ?: {groupElement(params, vec_E )}')
        print(f'Does vec_a have n ({n}) components in Z_q?: {len(vec_a) == n}')
        print(f'Are r, b, s and tau elements of Z_q ?: {0 <= r <= q and 0 <= b <= q and 0 <= s <= q and 0<=tau <= q }')
        print(f'Is c_B_m == com_ck(0;0)? : {com_B_k[m] == 1}')


        print(f'Is E_m == C? : {vec_E[m] == singleCiphertext}')

        # Now for equation verifications
        com_a_pow_x = bilinearMapCom(com_A,vec_x_prime, p)
        com_a0_times_a_pow_x = (com_A0 * com_a_pow_x)%p
        com_a_r = pc.commit(vec_a, params, r)
        print(f'Is com_a0 * com_a^x == com_a_r?: {com_a0_times_a_pow_x == com_a_r}')

        com_b_s = pc.commit(b, params, s)
        com_B_k_pow_x = bilinearMapCom(com_B_k,temp_vec_x, p) 
        print(f'Is com(b;s) == com_B_k^(x^k) ?: {com_B_k_pow_x== com_b_s}')

        e_k_power_x = bilinearMap(vec_E, temp_vec_x, p)
        enc_k_temp = elG._encrypt(pow(elG.publickey().g, b, p), tau)
        prod = [1,1]
        temp_vec_x_2 = [pow(x_multi, m-i-1,q) for i in range(m)]

        # Multiply each component of the power of x with each entry of vector a 
        result = [ (a * b)%q for a in temp_vec_x_2 for b in vec_a]
        ''' Result should contain a vector of length n*m. We need to split this 
            into components of n again to get m vectors in total
            Each of which is of size n'''
        pow_rhs = split_list(result, n)

        for i in range(m):
             tmp = bilinearMap(mat_c[i], pow_rhs[i], p)
             prod = [(x*y)%p for x,y in zip(prod, tmp)]
        # prod should contain C_i^(x^(m−i)*a).
        rhs_eq = [(i*j)%p for i,j in zip(enc_k_temp, prod)]
        print(f'Does the product of E_k`s raised to the power x^k == the right hand side ?: {e_k_power_x == rhs_eq}')
        end = timer()
        self.verifierTime.append(end-start)
        if (len(vec_a) == n and 0 <= r <= q and 0 <= b <= q and 0 <= s <= q and 0<=tau <= q and com_B_k[m] == 1 and vec_E[m] == singleCiphertext and com_a0_times_a_pow_x == com_a_r and com_B_k_pow_x== com_b_s and e_k_power_x == rhs_eq) == False:
            raise Exception("Error in Multi exponentiation argument verification")

    def multi_expo_optimization_arg_verify(self, pc, params, mew, singleCiphertext,  com_vec_b, vec_E, b, s, x_multi_optim ):
        start = timer()
        q, _, _, p = params 

        vec_x_prime = [pow(x_multi_optim, i, q) for i in range(2*mew - 1)]

        print(f'Does com_b have 2*mew - 1 {2*mew - 1} elements ?: {len(com_vec_b) == 2*mew - 1} ')
        print(f'Are all E_k`s elements of H: G x G ?: {groupElement(params, vec_E )}')

        print(f'Are b, s elements of Z_q? : (b: {0<=b<=q-1}) (s: {0<=s<=q-1}) ')

        print(f'Is c_B_{mew-1} == com_ck(0;0)? : {com_vec_b[mew-1] == 1}')

        print(f'Is E_{mew-1} == C? : {vec_E[mew-1] == singleCiphertext}')

        com_b_s = pc.commit(b, params, s)
        com_b_pow_x = bilinearMapCom(com_vec_b,vec_x_prime, p) 
        print(f'Is com(b;s) == com_b^(x) ?: {com_b_pow_x== com_b_s}')

        if (len(com_vec_b) == 2*mew - 1 and 0 <= b<= q-1 and 0 <= s <= q-1 and com_vec_b[mew-1] == 1 and vec_E[mew-1] == singleCiphertext and com_b_pow_x== com_b_s) == False:
            raise Exception("Error in Multi exponentiation optimization verification")

    def verifier_time(self):
        return sum(self.verifierTime)
class Prover:
    def __init__(self) -> None:
        self.proverTime = []
        self.permutation_init = 0
        self.rho_init = 0 

        self.vec_a_r1 = 0   # {pi(1), ... pi(N)}
        self.vec_r_r1 = 0       # Corresponding commitment randomizer for vec_a

        self.vec_b_r3 = 0   #{x^pi(i)}
        self.vec_s_r3 = 0   # Corresponding commitment randomizer for vec_b

        self.vec_d_r5 = 0   #{y*pi(i) + x^pi(i)}
        self.vec_t_r5 = 0   # Corresponding commitment randomizer for vec_d (& for vec_d_minus_z_r5 since opening for vec_z was vec_0)
        self.vec_d_minus_z_r5 = 0

        self.vec_b_prod_arg_r1 = 0  # { [prod(a_{1j}), prod(a_{2j}), ... , prod(a_{nj})}; product of each row of mat(a) = mat(d-z)
        self.s_prod_arg_r1 = 0      # Corresponding commitment randomizer for vec_b of the product argument

        self.mat_B_hadamard_prod_arg_r1 = 0     # [b_1, b_2, b_3, b_4] = [a_1, a_1a_2, ... , a_1a_2a_3a_4 ] where a_i are columns of vec_d_minus_z_r5 (i.e. a_1 = [d_1 - z, ..., d_4 - z])
        self.vec_s_hadamard_prod_arg_r1 = 0     # Corresponding commitment randomizer for mat_B of the Hadamard product argument

        self.mat_A_hadamard_prod_arg_r3 = 0     # [a_2, a_3, a_4, -1] i.e. columns 2 .. m of vec(d_i - z) with column of [-1, -1, -1, -1] as last entry 
        self.vec_r_hadamard_prod_arg_r3 = 0     # Corresponding commitment randomizer for mat_A of the Hadamard product argument. Contains [t[1], t[2], t[3], 0]
        self.mat_D_hadamard_prod_arg_r3 = 0     # [xb_1, x^2b_2, x^3b_3, xb_2 + x^2b_3 + x^3b_4]. i.e. first m-1 components are entries of mat_B of this argument * x^i and last component is vec_d
        self.vec_t_hadamard_prod_arg_r3 = 0     # Corresponding commitment randomizer for mat_D of the Hadamard product argument


        
        self.vec_a0_zero_arg_r1 = 0     # Vector of random values
        self.r0_zero_arg_r1 = 0         # Corresponding randomizer for vec_a0
        self.vec_bm_zero_arg_r1 = 0     # Vector of random values
        self.sm_zero_arg_r1 = 0         # Corresponding randomizer for vec_bm
        self.vec_r_zero_arg_r1 = 0      # [r0, vec_r_hadamard_prod_arg_r3]
        self.mat_A_zero_arg_r1 = 0      # [vec_a0, mat_A_hadamard_prod_arg_r3]
        self.mat_B_zero_arg_r1 = 0      # [mat_D_hadamard_prod_arg_r3, vec_bm]
        self.vec_s_zero_arg_r1 = 0      # [vec_t_hadamard_prod_arg_r3, vec_sm]
        self.vec_D_zero_arg_r1 = 0      # Contains the diagonal entries of the matrix given by A (*) B where A = mat_A_zero_arg_r1, B = mat_B_zero_arg_r1 and (*) is the bilinear map defined
        self.vec_t_zero_arg_r1 = 0      # Corresponding commitment randomizers used for vec_D of zero_arg. t[m+1] = 0 since corresponding diagonal entry d_{m+1} = 0


        self.vec_delta_single_val_arg_r1 = 0        # Vector of random values
        self.vec_b_single_val_arg_r1 = 0            # Vector of values [a_1, a_1a_2, a_1a_2a_3, a_1a_2a_3a_4 ] where a_i are entries of vec_b_prod_arg_r1
        self.vec_d_single_val_arg_r1 = 0            # Vector of random values
        self.r_d_single_val_arg_r1 = 0              # Random value used for commitment to vec_d_single_val_arg_r1
        self.vec_delta_d_single_val_arg_r1 = 0      # Vector used for com_delta
        self.s_1_single_val_arg_r1 = 0              # Random value used for commitment to vec_delta_d_single_val_arg_r1
        self.vec_delta_minus_single_val_arg_r1 = 0  # Vector used for com_triangle
        self.s_x_single_val_arg_r1 = 0              # Random value used for commitment to vec_delta_minus_single_val_arg_r1

        self.rho_multi_expo_arg_r0 = 0              # (-1) * rho_init *  vec_b_r3

        self.vec_a_0_multi_expo_arg_r1 = 0      # Vector of n random values
        self.r_0_multi_expo_arg_r1 = 0          # Corresponding randomizer for commitment to vec_a_0_multi_expo_arg_r1
        self.vec_b_multi_expo_arg_r1 = 0        # Vector of 2*m random values
        self.vec_s_multi_expo_arg_r1 = 0        # Corresponding vector of randomizers for commitment to vec_b_multi_expo_arg_r1
        self.vec_tau_multi_expo_arg_r1 = 0      # Vector of random values used for each ElGamal encryption

        self.vec_b_multi_expo_optim_r1 = 0      
        self.vec_s_multi_expo_optim_r1 = 0
        self.vec_tau_multi_expo_optim_r1 = 0

        self.mat_A_prime_multi_expo_optim_r3 = 0        # New secret values
        self.vec_r_prime_multi_expo_optim_r3 = 0        # Openings for ^
        self.rho_prime_multi_expo_optim_r3 = 0          # New rho 

        self.vec_a_0_multi_expo_mod_arg_r1 = 0
        self.r_0_multi_expo_mod_arg_r1 = 0
        self.vec_b_multi_expo_mod_arg_r1 = 0
        self.vec_s_multi_expo_mod_arg_r1 = 0
        self.vec_tau_multi_expo_mod_arg_r1 = 0

    def witness_initial(self):
        start = timer()
        
        permuted_list = random.sample(range(1, N+1), N)
        print(f'Permutation used: {permuted_list}')

        # Generating N random values for vector rho
        vec_rho = list_of_rand_values(N, q)
        print(f'Rho vector given by: {vec_rho}')
        self.permutation_init = permuted_list
        self.rho_init = vec_rho
        end = timer()
        self.proverTime.append(end-start)
        return permuted_list, vec_rho
    
    def round_1_computation(self, pc, params, m):
        start = timer()
        vec_r = list_of_rand_values(m, q)
        vec_a = self.permutation_init
        
        com_a = pc.generalCommit(params, vec_a, vec_r)

        self.vec_r_r1 = vec_r
        self.vec_a_r1 = vec_a
        end = timer()
        self.proverTime.append(end-start)
        return com_a

    def round_3_computation(self, pc, params, m, x_init):
        start = timer()
        q = params[0]
        vec_s = list_of_rand_values(m, q)
        vec_b = [pow(x_init, pi_i, q) for pi_i in self.vec_a_r1]
        # Uncomment below to make product arg fail (i.e. different permutation used)
        '''
        pi_mod = self.vec_a_r1
        pi_mod[-1] = 4
        vec_b = [pow(x_init, pi_i, q) for pi_i in pi_mod]
        '''
        com_b = pc.generalCommit(
            params, vec_b, vec_s)
        self.vec_s_r3 = vec_s
        self.vec_b_r3 = vec_b
        end = timer()
        self.proverTime.append(end-start)
        return com_b
    def round_5_computation(self, pc, params, m, n, y_init, vec_z):
        start = timer()
        vec_d_t1 = [(y_init * a_i) % q for a_i in self.permutation_init]
        vec_d = [(y_times_ai + b_i) % q for y_times_ai, b_i in zip(vec_d_t1, self.vec_b_r3)]
        print(f'Prover computes vector d  given by: y*a_i + b_i')
        
        
        vec_t_t1 = [(y_init * r_i) % q for r_i in self.vec_r_r1]
        vec_t = [(y_times_r + s_i) % q for y_times_r, s_i in zip(vec_t_t1, self.vec_s_r3)]
        
        print(f'Prover computes vector t given by: y*r_i + s_i')
        
        vec_d_minus_z = [(d_i+z_i) % q for d_i, z_i in zip(vec_d, vec_z)]
        com_d_minus_z = pc.generalCommit(params, vec_d_minus_z, vec_t)
        
        self.vec_d_r5 = vec_d
        self.vec_t_r5 = vec_t
        self.vec_d_minus_z_r5 = vec_d_minus_z
        end = timer()
        self.proverTime.append(end-start)
        return com_d_minus_z
        
    def product_arg_round_1_computation(self, pc, params):
        start = timer()
        q = params[0]
        s = secretsGenerator.randint(1, q-1)
        # s = random.randint(1, q-1)
        matA = split_list(self.vec_d_minus_z_r5, n)
        # Need to commit to product of each row of matA 
        vec_b = split_cols_and_multiply(matA, q)
        com_b = pc.commit(vec_b, params , s) 

        self.vec_b_prod_arg_r1 = vec_b
        self.s_prod_arg_r1 = s
        end = timer()
        self.proverTime.append(end-start)
        return com_b 

    def hadamard_prod_arg_round_1_computation(self, pc, params, m, n, q):
        start = timer()
        matA = split_list(self.vec_d_minus_z_r5, n) # Converts to a 2-d list of m vectors containing n values each 
        mat_B_i = []    # Shouuld contain all the vec_b_i's 
        vec_b_1 = matA[0]   # vec_a_1 essentially
        mat_B_i.append(vec_b_1)
        for i in range(1, len(matA)):
            col = matA[i]
            vec_b_i_temp = np.multiply(mat_B_i[i-1], col)  # Hadamard product; entry-wise multiplication
            mat_B_i.append( [(int(elt))%q for elt in vec_b_i_temp])   

        # Define s_i's 
        vec_s = list_of_rand_values(m, q)
        vec_s[0] = self.vec_t_r5[0]
        vec_s[-1] = self.s_prod_arg_r1

        vec_com_B = [pc.commit( mat_B_i[k], params, vec_s[k]) for k in range(m)]
        
        self.mat_B_hadamard_prod_arg_r1 = mat_B_i
        self.vec_s_hadamard_prod_arg_r1 = vec_s
        end = timer()
        self.proverTime.append(end-start)
        return vec_com_B
    
    def hadamard_prod_arg_round_3_computation(self, pc, params, m, n, q, x_had, vec_comA, vec_com_B):
        start = timer()
        matA = split_list(self.vec_d_minus_z_r5, n)

        temp_vec_x = [pow(x_had, i, q) for i in range(1, m)]

        # List of values from i = 1 .. m-1
        vec_com_D_i = [pow(vec_com_B[i], temp_vec_x[i], p) for i in range(len(vec_com_B) - 1)]

        #  cD =prod^{m−1}_{i=1} c^{xi}_{B_{i+1}}
        com_D_temp = [pow(vec_com_B[i+1], temp_vec_x[i], p) for i in range(len(vec_com_B) - 1)]

        com_D = reduce((lambda x, y: (x * y)%p), com_D_temp)
        vec_neg_1 = [-1] * n
        com_minus1 = pc.commit(vec_neg_1, params, 0 )
        # vectors a_2, a_3,...a_m
        mat_A_new = []
        mat_A_new.extend(matA[1:])
        mat_A_new.append(vec_neg_1)

        vec_comA_new = []
        vec_comA_new.extend(vec_comA[1:])
        vec_comA_new.append(com_minus1)

        vec_r_new = self.vec_t_r5[1:]
        vec_r_new.append(0)

        

        # vec_com_D_i already has m-1 elements 
        vec_com_D_new = []
        vec_com_D_new.extend(vec_com_D_i)
        vec_com_D_new.append(com_D)

        # Contains vec_d_i's up to m-1
        mat_D_i = []
        vec_d = []  # contains entries x^i * vec_b_{i+1}
        i = 0
        for col in self.mat_B_hadamard_prod_arg_r1:
            # Only need the first m-1 columns 
            if col != self.mat_B_hadamard_prod_arg_r1[0]:
                vec_d.append([(temp_vec_x[i-1] * elt) % q for elt in col])
                if col == self.mat_B_hadamard_prod_arg_r1[-1]:
                    break
            mat_D_i.append( [(temp_vec_x[i] * elt) % q for elt in col] )
            i += 1
        # mat_D should be {x*b_1 x^2*b_2 ... x^}
        # vec_d_final = x*b_2 + x^2*b_3 + ... x^{m-1}*b_m
        vec_d_final = [ (sum(x))%q for x in zip(*vec_d)]

        mat_D = []
        mat_D.extend(mat_D_i)
        mat_D.append(vec_d_final)

        # temp_vec_x already has m-1 elements (x^1, x^2, ..., x^{m-1})
        vec_t = np.multiply(temp_vec_x, self.vec_s_hadamard_prod_arg_r1[:-1]).tolist()

        t = sum(np.multiply(temp_vec_x, self.vec_s_hadamard_prod_arg_r1[1:]))
        vec_t.append(t)

        self.mat_A_hadamard_prod_arg_r3 = mat_A_new
        self.vec_r_hadamard_prod_arg_r3 = vec_r_new
        self.mat_D_hadamard_prod_arg_r3 = mat_D
        self.vec_t_hadamard_prod_arg_r3 = vec_t
        end = timer()
        self.proverTime.append(end-start)
        return vec_comA_new, vec_com_D_new

    def zero_arg_round_1_computation(self, pc, params, y_had):
        start = timer()
        q, g_i, h, p = params
        mat_A = []
        mat_A.extend(self.mat_A_hadamard_prod_arg_r3)
        vec_r = []
        vec_r.extend(self.vec_r_hadamard_prod_arg_r3)
        mat_B = []
        mat_B.extend(self.mat_D_hadamard_prod_arg_r3)
        vec_s = []
        vec_s.extend(self.vec_t_hadamard_prod_arg_r3)
 
        n = len(mat_A[0])
        m = len(mat_A)
        vec_a0 = list_of_rand_values(n, q)
        vec_bm = list_of_rand_values(n, q)
        r_0 = secretsGenerator.randint(1, q-1)
        s_m = secretsGenerator.randint(1, q-1)
        com_a0 = pc.commit(vec_a0, params, r_0)
        com_bm = pc.commit(vec_bm, params, s_m)
        ''' Create a new matrix with vec_a0 as the first entry follow by all other vectors of mat_A'''
        vec_r.insert(0, r_0)
        mat_A.insert(0, vec_a0)
        mat_B.insert(len(mat_B), vec_bm)
        vec_s.insert(len(vec_s), s_m)
        print(f'Dimensions of matrix A should be n x m + 1: {len(mat_A) == m+1}')

        vec_D = []
        # We need values d_0,...d_{2m} which are the diagonals of 
        for k in range(2*m + 1):
            # vec_d_k is used to keep track of which indices were used to compute the corresponding d_k
            vec_d_k = []
            d_k = 0
            for i in range(m+1):
                j = m-k + i 
                if 0 <= j <= m:
                    vec_d_k.append((i,j))
                    d_k = (d_k + zeroArgBilinearMap(mat_A[i], mat_B[j], y_had))%q
            vec_D.append((vec_d_k, d_k))

        print(f'vec_D should contain a total of 2m+1 elements ({2*m + 1}) : {len(vec_D)}')
        vec_t = list_of_rand_values(2*m+1, q)
        vec_t[m+1] = 0
        vec_com_D = [pc.commit(vec_D[i][1], params, vec_t[i]) for i in range(2*m+1)]

        self.vec_r_zero_arg_r1 = vec_r
        self.mat_A_zero_arg_r1 = mat_A
        self.mat_B_zero_arg_r1 = mat_B
        self.vec_s_zero_arg_r1 = vec_s
        self.vec_a0_zero_arg_r1 = vec_a0      
        self.vec_bm_zero_arg_r1 = vec_bm
        self.r0_zero_arg_r1 = r_0      
        self.sm_zero_arg_r1 = s_m
        self.vec_D_zero_arg_r1 = vec_D
        self.vec_t_zero_arg_r1 = vec_t
        end = timer()
        self.proverTime.append(end-start)
        return com_a0, com_bm, vec_com_D

    def zero_arg_round_3_computation(self, params, m, x_zeroArg, vec_x_m, vec_x_m_j, vec_x_m_k):
        start = timer()
        q, g_i, h, p = params
        vec_a = []
        for i in range(m+1):
            vec_a.append( [ (vec_x_m[i] * elt)%q for elt in self.mat_A_zero_arg_r1[i]])
        vec_a = [sum(x) for x in zip(*vec_a)]
        r = int(np.dot(vec_x_m, self.vec_r_zero_arg_r1))
        vec_b = []
        for i in range(m+1):
            vec_b.append( [ (vec_x_m_j[i] * elt)%q for elt in self.mat_B_zero_arg_r1[i]])
        vec_b = [sum(x) for x in zip(*vec_b)]

        s = int(np.dot(vec_x_m_j, self.vec_s_zero_arg_r1))
        t = int(np.dot(vec_x_m_k, self.vec_t_zero_arg_r1))
        end = timer()
        self.proverTime.append(end-start)
        return vec_a, vec_b, r, s, t 

    def single_val_arg_round_1_computation(self, pc, params):
        start = timer()
        vec_a = self.vec_b_prod_arg_r1 
        n = len(vec_a)
        q, g_i, h, p = params

        vec_b = []
        vec_b.append(vec_a[0])
        for i in range(1, len(vec_a)):
            vec_b.append((vec_b[i-1]*vec_a[i])%q)

        vec_d = list_of_rand_values(n, q)
        r_d = secretsGenerator.randint(1, q-1)
        vec_delta = list_of_rand_values(n, q)
        vec_delta[0] = vec_d[0]
        vec_delta[len(vec_delta)-1] = 0 

        s_1 = secretsGenerator.randint(1, q-1)
        s_x = secretsGenerator.randint(1, q-1)

        com_d = pc.commit(vec_d, params, r_d)

        # -delta_1*d_2, ... , -delta_{n-1} * d_n
        vec_delta_d = [ (delta_i * d_j * (-1))%q for delta_i, d_j in zip(vec_delta[:len(vec_delta) -1], vec_d[1:])]
        com_delta = pc.commit(vec_delta_d, params, s_1)

        vec_a_delta = [(a_i *  delta_j)%q for a_i, delta_j in zip(vec_a[1:], vec_delta[:len(vec_delta) -1])]
        vec_b_d = [(b_i *  d_j)%q for b_i, d_j in zip(vec_b[:len(vec_b)-1], vec_d[1:])]
        vec_delta_minus = [ (vec_delta_i - vec_a_delta_i - vec_b_d_i)%q for vec_delta_i, vec_a_delta_i, vec_b_d_i in zip(vec_delta[1:], vec_a_delta, vec_b_d) ]
        com_triangle = pc.commit(vec_delta_minus, params, s_x)

        self.vec_delta_single_val_arg_r1 = vec_delta
        self.vec_b_single_val_arg_r1 = vec_b
        self.vec_d_single_val_arg_r1 = vec_d 
        self.r_d_single_val_arg_r1 = r_d 
        self.vec_delta_d_single_val_arg_r1 = vec_delta_d
        self.s_1_single_val_arg_r1 = s_1
        self.vec_delta_minus_single_val_arg_r1 = vec_delta_minus
        self.s_x_single_val_arg_r1 = s_x
        end = timer()
        self.proverTime.append(end-start)
        return com_d, com_delta, com_triangle
    
    def single_val_arg_round_3_computation(self, pc, params, x_5_3):
        start = timer()
        q, g_i, h, p = params
        vec_a_bar = [ (x_5_3 * a_i + d_i)%q for a_i, d_i in zip(self.vec_b_prod_arg_r1, self.vec_d_single_val_arg_r1)]
        r_bar = (x_5_3*self.s_prod_arg_r1 + self.r_d_single_val_arg_r1)%q 
        vec_b_bar = [ (x_5_3 * b_i + delta_i)%q for b_i, delta_i in zip(self.vec_b_single_val_arg_r1, self.vec_delta_single_val_arg_r1)]
        s_bar = (x_5_3*self.s_x_single_val_arg_r1 + self.s_1_single_val_arg_r1)%q
        end = timer()
        self.proverTime.append(end-start)
        return vec_a_bar, vec_b_bar, r_bar, s_bar
    def multi_exponentiation_round_0_computation(self, params, x_init, N):
        start = timer()
        q, _, _, p = params 

        # rho = - vec_rho * vec_b
        rho_tmp = (sum([(rho_i*b_i) % q for rho_i, b_i in zip(self.rho_init, self.vec_b_r3) ])) % q
        rho = (rho_tmp*(-1))%q
        vec_x = [pow(x_init, i, q) for i in range(1, N+1)]

        self.rho_multi_expo_arg_r0 = rho 
        end = timer()
        self.proverTime.append(end-start)
        return rho, vec_x

    def multi_exponentiation_round_1_computation(self, pc, elG, params, m, n, shuffled_ciphertexts, rho):
        start = timer()
        q, g_i, h, p= params
        mat_a_j = split_list(self.vec_b_r3, n)

        # Split the ciphertexts into m vectors of length n each  
        mat_c = split_list(shuffled_ciphertexts, n )

        vec_a_0 = list_of_rand_values(n, q)
        r_0 = secretsGenerator.randint(1, q-1)
        # Compute b_0, s_0, tau_0, ... , b_2m-1, s_2m-1, tau_2m-1 from Z_q
        vec_b = list_of_rand_values(2*m, q)
        vec_s = list_of_rand_values(2*m, q)
        vec_tau = list_of_rand_values(2*m, q)
        # set bm =0, sm = 0, τm = ρ
        vec_b[m] = 0
        vec_s[m] = 0
        vec_tau[m] = rho 
        com_A0 = pc.commit(vec_a_0, params, r_0)
        vec_com_B_k = pc.generalCommit(params, vec_b, vec_s)
        

        mat_A = [vec_a_0] + mat_a_j
        vec_E = []
        diagonal_vals = []
        print(f'Calculating E_k`s')
        # For each E_k compute the ElGamal encryption part  
        for k in range(2*m):
            enc_k = elG._encrypt(pow(elG.publickey().g, vec_b[k], p), vec_tau[k])
            # 1 <= i <= m
            vec_E_k = []
            E_k = [1,1] 
            for i in range(0,m):
                j = k-m + i+1       # Plus one because of index issue
                if 0 <= j <=m:
                    c_i_pow_a = bilinearMap(mat_c[i], mat_A[j], p)
                    vec_E_k.append(c_i_pow_a)
                    E_k =  [(x*y)%p for x,y in zip(E_k,c_i_pow_a)]
            '''enc_k is the encryption part. vec_E_k contains all the values that will be part of the diagonal
                E_k should contain the product of all the values in vec_E_k
                vec_E should contain the final result E_k for k = 0 ... 2m-1. All values before vec_E are just for 
                visual purposes'''
            diagonal_vals.append(vec_E_k)
            vec_E.append(list(tuple(
            (i*j) % p for i, j in zip(E_k, enc_k))))
        print(f'Finished calculating E_k`s')
        self.vec_a_0_multi_expo_arg_r1 = vec_a_0
        self.r_0_multi_expo_arg_r1 = r_0
        self.vec_b_multi_expo_arg_r1 = vec_b
        self.vec_s_multi_expo_arg_r1 = vec_s
        self.vec_tau_multi_expo_arg_r1 = vec_tau
        
        end = timer()
        self.proverTime.append(end-start)
        return com_A0, vec_com_B_k, vec_E

    def multi_exponentiation_round_3_computation(self, params, m, n, x_multi):
        start = timer()
        q, _, _, p = params 
        vec_x_prime = [pow(x_multi, i, q) for i in range(1, m+1)]

        vec_a_j = split_list(self.vec_b_r3, n)      # Feel like this computation is done too often. Chack which elt corresponds to it 

        vec_Ax = []
        for i in range(len(vec_a_j[0])):
            resultTemp = 0
            for j in range(len(vec_a_j)):
                resultTemp = ( resultTemp + ((vec_a_j[j][i]*vec_x_prime[j])%q) )%q
            vec_Ax.append(resultTemp)
        
        vec_a = [(a+b)%q for a,b in zip(self.vec_a_0_multi_expo_arg_r1, vec_Ax)]

        # r = r0+~r·x
        '''Type casting from numpy object back to integers'''
        r = int((self.r_0_multi_expo_arg_r1 + np.dot(self.vec_s_r3, vec_x_prime))%q)

        temp_vec_x = [pow(x_multi, i, q) for i in range(2*m)]
        b = int(np.dot(self.vec_b_multi_expo_arg_r1, temp_vec_x)%q)
        s = int(np.dot(self.vec_s_multi_expo_arg_r1, temp_vec_x)%q)
        tau = int(np.dot(self.vec_tau_multi_expo_arg_r1, temp_vec_x)%q)
        end = timer()
        self.proverTime.append(end-start)
        return vec_a, r, b, s, tau 

    def multi_expo_optimization_round_1_computation(self, pc, elG, params, m_prime, n, shuffled_ciphertexts, rho, mew):
        start = timer()
        q, g_i, h, p= params
        mat_A = split_list(self.vec_b_r3, n)            # original challenge x raised to permutation used

        # Split the ciphertexts into m vectors of length n each  
        mat_c = split_list(shuffled_ciphertexts, n )

        # Compute b_0, s_0, tau_0, ... , b_2mew-2, s_2mew-2, tau_2mew-2 from Z_q
        vec_b = list_of_rand_values(2*mew - 1, q)
        vec_s = list_of_rand_values(2*mew -1, q)
        vec_tau = list_of_rand_values(2*mew -1 , q)

        # set b_{mew-1} =0, s_{mew-1} = 0, τ_{mew-1} = ρ
        vec_b[mew - 1] = 0
        vec_s[mew - 1] = 0
        vec_tau[mew - 1] = rho 
       
        vec_com_B_k = pc.generalCommit(params, vec_b, vec_s)
        vec_E = []
        diagonal_vals = [] 
        print(f'Calculating E_k`s')  
        for k in range(2*mew - 1):
            enc_k = elG._encrypt(pow(elG.publickey().g, vec_b[k], p), vec_tau[k])
            vec_E_k = []
            E_k = [1,1] 
            for l in range(0, m_prime):
                for i in range(1,mew+1):
                    j = k+1-mew +i 
                    if 1 <= j <=mew:
                        c_i_pow_a = bilinearMap(mat_c[mew * l + i - 1], mat_A[mew*l+j-1], p)        # Should be the product of n elements raised to each power (output: a single ciphertext)
                        vec_E_k.append((mew * l + i, mew*l+j ))
                        E_k =  [(x*y)%p for x,y in zip(E_k,c_i_pow_a)]
                '''enc_k is the encryption part. vec_E_k contains the indices i,j for C_i^a^j that will be part of the diagonal
                    E_k should contain the product of all the values in vec_E_k
                    vec_E should contain the final result E_k for k = 0 ... 2m-1. All values before vec_E are just for 
                    visual purposes'''
            diagonal_vals.append(vec_E_k)
            vec_E.append(list(tuple((i*j) % p for i, j in zip(E_k, enc_k))))
    
        self.vec_b_multi_expo_optim_r1 = vec_b
        self.vec_s_multi_expo_optim_r1 = vec_s
        self.vec_tau_multi_expo_optim_r1 = vec_tau
        
        end = timer()
        self.proverTime.append(end-start)
        return vec_com_B_k, vec_E
    
    def multi_expo_optimization_round_3_computation(self, params, m_prime, n, x_multi_optim, mew, shuffled_ciphertexts, vec_comA, vec_E):
        start = timer()
        q, _, _, p = params 
        vec_x_prime = [pow(x_multi_optim, i, q) for i in range(0, 2*mew-1)]

        b = int(np.dot(self.vec_b_multi_expo_optim_r1, vec_x_prime)%q)
        s = int(np.dot(self.vec_s_multi_expo_optim_r1, vec_x_prime)%q)
        rho_prime = int(np.dot(self.vec_tau_multi_expo_optim_r1, vec_x_prime)%q)

        mat_c = split_list(shuffled_ciphertexts, n )

        mat_A = split_list(self.vec_b_r3, n)      # Feel like this computation is done too often. Chack which elt corresponds to it 

        mat_A_prime = []    # Each entry is a vector 
        vec_r_prime = []
        entries_used = []
        entries_used_for_new_ciph = []
        new_mat_c = []      # New vectors of ciphertexts
        com_A_prime = []
        for l in range(1, m_prime + 1):
            resultTemp = [0]*len(mat_A[0])
            vec_r_tmp = 0
            entries_used_tmp = []
            entries_used_tmp2 = []
            result_ciph = [(1,1)]*len(mat_c[0]) 
            com_result = 1
            for j in range(1, mew+1):
                index = mew*(l-1) + j -1 # -1 because a_i is at ith index
                tmp_vector = [ (vec_x_prime[j-1] * val)%q for val in mat_A[index] ]
                resultTemp = [(a+b)%q for a,b in zip(resultTemp, tmp_vector)]
                entries_used_tmp.append((j-1, index+1))
                vec_r_prime_component = (vec_x_prime[j-1] * self.vec_s_r3[index])%q 
                vec_r_tmp = (vec_r_tmp + vec_r_prime_component)%q 
                # For new vectors of ciphertexts and commitments
                entries_used_tmp2.append((mew - j, index+1))
                x_component_ciph = [vec_x_prime[mew-j]] * n      # Multiply each element of mat_C[i] which contains n elements with the same power
                result_tmp = ciphertextPower(mat_c[index], x_component_ciph, p)
                result_ciph = [ ((a*c)%p,(b*d)%p) for (a,b), (c,d) in zip(result_ciph, result_tmp)]       # Cool
                com_result = (com_result * pow(vec_comA[index], vec_x_prime[j-1], p))%p
            
            entries_used_for_new_ciph.append(entries_used_tmp2)
            new_mat_c.append(result_ciph)
            com_A_prime.append(com_result)
            entries_used.append(entries_used_tmp)
            mat_A_prime.append(resultTemp)
            vec_r_prime.append(vec_r_tmp)

        exponent = (b*-1)+q
        enc = elG._encrypt(pow(elG.publickey().g, exponent, p), 0)
        E_power_x = bilinearMap(vec_E, vec_x_prime, p)
        single_ciphertext_prime = list(tuple((i*j) % p for i, j in zip(E_power_x, enc)))

        self.mat_A_prime_multi_expo_optim_r3 = mat_A_prime
        self.vec_r_prime_multi_expo_optim_r3 = vec_r_prime
        self.rho_prime_multi_expo_optim_r3 = rho_prime

        end = timer()
        
        self.proverTime.append(end-start)
        return b,s, single_ciphertext_prime, new_mat_c, com_A_prime
    
    def multi_exponentiation_modified_round_1_computation(self, pc, elG, params, m, n, shuffled_ciphertexts, rho):
        start = timer()
        q, g_i, h, p= params
        mat_a_j = self.mat_A_prime_multi_expo_optim_r3
        mat_c = shuffled_ciphertexts   # shuffled_ciphertexts already a 2d matrix

        vec_a_0 = list_of_rand_values(n, q)
        r_0 = secretsGenerator.randint(1, q-1)
        vec_b = list_of_rand_values(2*m, q)
        vec_s = list_of_rand_values(2*m, q)
        vec_tau = list_of_rand_values(2*m, q)
        # set bm =0, sm = 0, τm = ρ
        vec_b[m] = 0
        vec_s[m] = 0
        vec_tau[m] = rho 
        com_A0 = pc.commit(vec_a_0, params, r_0)
        vec_com_B_k = pc.generalCommit(params, vec_b, vec_s)
        
        mat_A = [vec_a_0] + mat_a_j
        vec_E = []
        diagonal_vals = []
        print(f'Calculating E_k`s')
        # For each E_k compute the ElGamal encryption part  
        for k in range(2*m):
            enc_k = elG._encrypt(pow(elG.publickey().g, vec_b[k], p), vec_tau[k])
            # 1 <= i <= m
            # Contains the corresponding elements of E_k
            vec_E_k = []
            E_k = [1,1] 
            for i in range(0,m):
                j = k-m + i+1       # Plus one because of index issue 
                # Only when this condition is satisfied should we add C_i^a_j
                if 0 <= j <=m:
                    c_i_pow_a = bilinearMap(mat_c[i], mat_A[j], p)
                    vec_E_k.append((i+1, j))
                    E_k =  [(x*y)%p for x,y in zip(E_k,c_i_pow_a)]
            '''enc_k is the encryption part. vec_E_k contains all the values that will be part of the diagonal
                E_k should contain the product of all the values in vec_E_k
                vec_E should contain the final result E_k for k = 0 ... 2m-1. All values before vec_E are just for 
                visual purposes'''
            diagonal_vals.append(vec_E_k)
            vec_E.append(list(tuple(
            (i*j) % p for i, j in zip(E_k, enc_k))))

        self.vec_a_0_multi_expo_mod_arg_r1 = vec_a_0
        self.r_0_multi_expo_mod_arg_r1 = r_0
        self.vec_b_multi_expo_mod_arg_r1 = vec_b
        self.vec_s_multi_expo_mod_arg_r1 = vec_s
        self.vec_tau_multi_expo_mod_arg_r1 = vec_tau

        end = timer()
        self.proverTime.append(end-start)
        return com_A0, vec_com_B_k, vec_E
    
    def multi_exponentiation_modified_round_3_computation(self, params, m_prime, n, x_multi):
        start = timer()
        q, _, _, p = params 
        vec_x_prime = [pow(x_multi, i, q) for i in range(1, m_prime+1)]
        vec_a_j = self.mat_A_prime_multi_expo_optim_r3
        vec_Ax = []
        for i in range(len(vec_a_j[0])):
            resultTemp = 0
            for j in range(len(vec_a_j)):
                resultTemp = ( resultTemp + ((vec_a_j[j][i]*vec_x_prime[j])%q) )%q
            vec_Ax.append(resultTemp)
        
        vec_a = [(a+b)%q for a,b in zip(self.vec_a_0_multi_expo_mod_arg_r1, vec_Ax)]
        ''' r = r0+~r·x
            Type casting from numpy object back to integers'''
        r = int((self.r_0_multi_expo_mod_arg_r1 + np.dot(self.vec_r_prime_multi_expo_optim_r3, vec_x_prime))%q)


        temp_vec_x = [pow(x_multi, i, q) for i in range(2*m_prime)]
        b = int(np.dot(self.vec_b_multi_expo_mod_arg_r1, temp_vec_x)%q)
        s = int(np.dot(self.vec_s_multi_expo_mod_arg_r1, temp_vec_x)%q)
        tau = int(np.dot(self.vec_tau_multi_expo_mod_arg_r1, temp_vec_x)%q)
        end = timer()
        self.proverTime.append(end-start)
        return vec_a, r, b, s, tau 
    
    def get_rho(self):
        return self.rho_prime_multi_expo_optim_r3
    def prover_time(self):
        return sum(self.proverTime)
    
'''The verifier (or third party) generates the public parameters for the commitment scheme
    We'll just stick to a third party generating it '''
class PedersenCommitment:
    
    def setup(self, noOfBits_modulus, noOfBits_subgroup, n):

        q = secrets.randbits(noOfBits_subgroup)
        while(number.isPrime(q) == False):
            q = secrets.randbits(noOfBits_subgroup)

    
        print(f'Prime q: {q}')
        print(f'Length of q (in bits): {q.bit_length()}')
        r = pow(2, noOfBits_modulus-noOfBits_subgroup)
        while True:
            p = r*q + 1
            if number.isPrime(p):
                print("p = ", p)
                print(f'Length of p (in bits): {p.bit_length()}')
                break
            r += 1
        
        g_i = []
        
        '''
        The verifier can find a generator g of G by computing g = g_1^{\frac{p-1}{q}}mod p, where g_1 is a generator of the group Z_p.
        We know what the group Z_p is but we don't know G (which is a subgroup of Z_p* of order q)
        Z_p = {0, 1, ..., p-1}

        '''
        g_1 = secretsGenerator.randint(2, p-1)
        while pow(g_1, r, p) == 1:
             g_1 = secretsGenerator.randint(2, p-1)
        
        g = pow(g_1, r, p)
        g_i.append(g)

        randomPowers = []
        while len(randomPowers) != n:
            x = secretsGenerator.randint(2, q-1)
            if x not in randomPowers and x != g_1:
                randomPowers.append(x)

        for val in randomPowers:
            g_i.append(pow(g, val, p))
        h = g_i[-1] # Since g_i contains n+1 elts in total

        if len(g_i) != len(set(g_i)):
            raise Exception("Values of the generators are not unique")
        
        g_i = g_i[:-1]  # To get rid of h

        print(f'Length of g_i should be equal to n ({n}): {len(g_i)}')
        
        return q, g_i, h, p
    
    # Commit to a vector of m <= n values and randomness r
    def commit(self, m, params, r):
        q, g_i, h, p = params
        prod_g = 1
        # Should take care of the case for com_ck(0;0)
        if type(m)!=list:
            m = [m]
        for i in range(len(m)):
            prod_g = (prod_g * pow(g_i[i], m[i], p))%p
        c = (pow(h, r, p) * prod_g) % p
        return c

    '''
        vec_a: dimension N (couuld also be a (n x m) matrix with m columns each of size n such that N = nm)
        vec_r: dimension m
        Output: Return an element of G^m such that:
            com_ck(~a; ~r) = (com_ck(a_1, ..., a_n; r_1), ... ,  com_ck(a_{(m-1)n+1}, ..., a_N; r_m)
    '''
    def generalCommit(self, params, vec_a, vec_r):
        # Compute the value n
        n = len(vec_a) // len(vec_r)
        # Compute each column of the matrix representation of vec_a
        columns = split_list(vec_a, n)
        resultantCommitment = []
        for elt in range(len(columns)):
            resultantCommitment.append(
                self.commit(columns[elt], params, vec_r[elt]))
        return resultantCommitment

    def open(self, param, c, m, *r):
        q, g_i, h, p = param

        rSum = 0
        for rEl in r:
            rSum += rEl
        prod_pow = 1
        for msg in range(len(m)):
            prod_pow = (prod_pow * pow(g_i[msg], m[msg], p)) % p
        return c == (prod_pow * pow(h, rSum, p)) % p


    def add(self, param, *c):
        p = param[3]

        cSum = 1
        for cEl in c:
            cSum *= cEl
        return cSum % p

    def genAdd(self, param, c1, c2):
        p = param[3]
        genResult = []
        for i, j in zip(c1, c2):
            genResult.append(self.add(param, i, j))
        return genResult


def split_list(lst, n):
    return [lst[i:i + n] for i in range(0, len(lst), n)]

# Since each element is a generator and the order of G is q
def groupElement(params, ciphertexts):
    q, _, _, p = params
    for ciphertext in ciphertexts:
        if pow(ciphertext[0],q, p) != 1 or pow(ciphertext[1], q, p) != 1:
            raise Exception(f"Ciphertext {(ciphertext[0], ciphertext[1])} is not a group element")
    return True 

def list_of_rand_values(n, q):
    result = []
    while len(result) != n:
        rndVal = secretsGenerator.randint(1, q-1)
        if rndVal not in result:
            result.append(rndVal)
    return result

def isComGroupElement(params, com):
    q, _, _, p = params
    if pow(com,q,p) != 1:
        raise Exception(f"Commitment {com} is not a group element")
    return True

def isGenComGroupElement(params, com):
    q, _, _, p = params
    for val in com:
        if pow(val, q, p) != 1:
            raise Exception(f'Value: {val} not a group element')
    return True

# Raising each component of a vector of ciphertexts to the same power
# Output: n length ciphertext
def ciphertextPower(cipherText, vec_x, p):
    result = []
    for i in range(len(cipherText)):
        result.append([pow(cipherText[i][0], vec_x[i], p), pow(
            cipherText[i][1], vec_x[i], p)])
    return result 


'''TODO: Can probably make this more efficient '''
# Raising a vector of ciphertexts ( each a tuple) to a power mod p
def bilinearMap(cipherText, vec_x, p):
    result = []
    for i in range(len(cipherText)):
        result.append([pow(cipherText[i][0], vec_x[i], p), pow(
            cipherText[i][1], vec_x[i], p)])
    finalResulta = 1
    finalResultb = 1
    for j in range(len(result)):
        finalResulta = (finalResulta * result[j][0]) % p
        finalResultb = (finalResultb * result[j][1]) % p

    return [finalResulta, finalResultb]

# Map for commitment raised to a power. Takes care of a single commitment as well as a vector of commitments
def bilinearMapCom(commitment, vec_x, p):
    result = 1
    if type(commitment) != list:
        commitment = [commitment]
    for i in range(len(commitment)):
        result = (result * (pow(commitment[i], vec_x[i], p))) %p
    return result 

# Given a 2-d list, where each sub-list corresponds to a column, multiply all entries of each row together
def split_cols_and_multiply(matrix, q):
    result = []
    for col in range(len(matrix[0])):
        row_vec = [row[col] for row in matrix]
        rowProduct = reduce(lambda x, y: (x*y) % q, row_vec)
        result.append(rowProduct)
    return result 

'''
    Given ciphertexts C_11, . . . , C_mn and a single ciphertext C, we will give in this section an argument of knowledge of openings 
    of commitments c_A to A = {aij} such that C = E_pk(1; rho)C_i^{a_i} and c_A = com_ck(A; r) 
    where C_i = (C{i1}, . . . , C{in}) and a_j = (a_{1j} , . . . , a_{nj})

    So a matrix of the ciphertexts would look like:
    [C_11  C_12  C_13  ...  C_1n]           [C_11  C_12  C_13  C_14]        [C_1   C_2   C_3   C_4  ]
    | .     .     .    ...   .  |           [C_21  C_22  C_23  C_24]        [C_5   C_6   C_7   C_8  ]
    | .     .     .    ...   .  |       =   [C_31  C_22  C_33  C_34]   =    [C_9   C_10  C_11  C_12 ]
    [C_n1  C_n2  C_n3  ...  C_nm]           [C_41  C_22  C_43  C_44]        [C_13  C_14  C_15  C_16 ]

    Given shuffled_ciphertexts = {C_1, C_2, C_3, ..., C_14, C_15, C_16}:
        Split_list will create a 2-D list mat_C: 
           [ [C_1   C_2   C_3   C_4  ] 
             [C_5   C_6   C_7   C_8  ]
             [C_9   C_10  C_11  C_12 ]
             [C_13  C_14  C_15  C_16 ] ]
            
    '''
def multiExpoArgument(pc, elG, m, n, params, prover, verifier, org_ciphertexts, shuffled_ciphertexts, com_A, x_init ):
    start = timer()
    q, _, _, p = params 
    N = m*n
    rho, vec_x = prover.multi_exponentiation_round_0_computation(params, x_init, N)

    singleCiphertext = bilinearMap(org_ciphertexts, vec_x, p)
    print(f'Ciphertexts raised to vector x: {singleCiphertext}')

    print(f'Multi-Exponentiation Argument:')

    print(f'Round 1: Initial Message')
    com_A0, com_B_k, vec_E = prover.multi_exponentiation_round_1_computation(pc, elG, params, m, n, shuffled_ciphertexts, rho)

    print(f'Round 2: Verifier sends a challenge x')
    x_multi = verifier.generateChallenge(params)
    print(f'Challenge sent by verifier for MEA : {x_multi}')

    print(f'Round 3: Prover calculates and sends the respective openings')
    vec_a, r, b, s, tau  = prover.multi_exponentiation_round_3_computation(params, m, n, x_multi)

    print(f'Round 4: Verification')
    verifier.multi_expo_arg_verify(pc, params, m, n, com_A0, com_B_k, vec_E, vec_a, r, b, s, tau, x_multi, singleCiphertext, com_A, shuffled_ciphertexts )
    print(f'End of Multi Exponentiation Argument')
    end = timer()
    print(f'Total time taken for multi-exponentiation argument: {end - start}')

def multiExpoArgumentOptimized(pc, elG, m_prime, n, params, prover, verifier, singleCiphertext, shuffled_ciphertexts, com_A_prime, mew ):
    start = timer()
    rho = prover.get_rho()      # rho' (rho_prime)

    print(f'Multi-Exponentiation (optimized) Argument:')
    print(f'Round 1: Initial Message')
    com_A0, com_B_k, vec_E = prover.multi_exponentiation_modified_round_1_computation( pc, elG, params, m_prime, n, shuffled_ciphertexts, rho)
    
    print(f'Round 2: Verifier sends a challenge x')
    x_multi = verifier.generateChallenge(params)
    print(f'Challenge sent by verifier for MEA optimized: {x_multi}')

    print(f'Round 3: Prover calculates and sends the respective openings')
    vec_a, r, b, s, tau  = prover.multi_exponentiation_modified_round_3_computation(params, m_prime, n, x_multi)

    print(f'Round 4: Verification')
    verifier.multi_expo_optim_arg_verify(pc, params, m_prime, n, com_A0, com_B_k, vec_E, vec_a, r, b, s, tau, x_multi, singleCiphertext, com_A_prime, shuffled_ciphertexts )
    print(f'End of Multi Exponentiation Argument')
    end = timer()
    print(f'Total time taken for multi-exponentiation argument: {end - start}')

# bilinear map ∗ : Z^n_q × Z^n_q → Z_q by (a_1, . . . , a_n)T ∗ (d_1, . . . , d_n)T =sum^n_{j=1} a_{j} \times d_{j} \times y^{j}
def zeroArgBilinearMap(vec_a, vec_d, y):
    return sum([ (vec_a[i]*vec_d[i]*pow(y, i+1, q))%q for i in range(len(vec_a))])
     
'''
    b = (yi + x^i - z) should be the product of all the vals in the matrix A (d_i - z)
    A is a matrix of the values d_i - z (where d_i := y*pi(i) + x^{pi(i)}) (Elements of Z_q)
    vec_r is the vector used for generating comA
    vec_comA is an element of G^m
    b is an element of Z_q
    Prover's witness: A (Z^{n x m}), vec_r (Z_q^m)
'''
def productArgument(pc, elG, m, n, params, prover, verifier, b, vec_comA):
    start = timer()
    print(f'Round 1: Initial message')
    com_b = prover.product_arg_round_1_computation(pc, params)
    # Engage in an SHVZK argument of knowledge as described in Section 5.1
    hadamardProdArgument(pc, elG, m, params, prover, verifier, vec_comA, com_b)
    singleValueProdArgument(pc, params, com_b, b)
    print(f'Round 2: Verification ')
    verifier.product_arg_verify(params, com_b)
    print(f'End of product argument')
    end = timer()
    print(f'Total time taken for Product argument: {end - start}')

'''
com_b contains commitment to product of values of rows of the matrix below
        Row 1 contains values [ d_1 - z, d_5 - z, d_9 - z, d_13 -z ]
                                .
                                .
        Row 4 contains values [ d_4 - z, d_8 - z, d_12 - z, d_16 -z ]
    Each entry of vec_b is a single value (product of row values)
vec_comA is the commitment to the matrix containing values d_i - z    
'''
def hadamardProdArgument(pc, elG, m, params, prover, verifier, vec_comA, com_b):
    q, g_i, h, p = params
    
    print(f'Engaging in Hadamard Product Argument:')
    print(f'Round 1: Prover performs the initial computations')
    vec_com_B = prover.hadamard_prod_arg_round_1_computation(pc, params, m, n, q)

    print(f'Round 2: Verifier sends two challenges: ', end=" ")

    # Challenge 
    x_had = verifier.generateChallenge(params)
    y_had = verifier.generateChallenge(params)
    print(f'x_had: {x_had}, y_had: {y_had}')

    print(f'Round 3: Prover performs computations and engages in a zero argument for the committed values')    
    vec_comA_new, vec_com_D_new = prover.hadamard_prod_arg_round_3_computation(pc, params, m, n, q, x_had, vec_comA, vec_com_B)

    print(f'Engage in a zero argument for the committed values defined during this round')
    zeroArgument(pc, elG, params, y_had, vec_comA_new, vec_com_D_new)

    print(f'Round 4: Verification')
    verifier.hadamard_arg_verify(params, vec_comA, com_b, vec_com_B )
    print(f'Hadamard product argument finished')

'''
    mat_A = [a_2, ... , a_m, -1] = [(d_5 - z)   (d_9 - z)   (d_13 - z)   (-1)]
        Corresponds to: mat_A_hadamard_prod_arg_r3
    vec_r = Randomizers used for mat_A
        Corresponds to: vec_r_hadamard_prod_arg_r3
    mat_B = [b_1 b_2 ... b_m] = [d_1 d_2 .. d_{m-1} d] = [ x(a_1) x^2(a_1a_2) ....  x^{m-1}(a_1a_2...a_m-1) (x(a_1) + x^2(a_1)(a_2) + ... )]
        Corresponds to: mat_D_hadamard_prod_arg_r3
    vec_s = Randomizers used for mat_B
        Corresponds to: vec_t_hadamard_prod_arg_r3
    '''
def zeroArgument(pc, elG, params, y_had, vec_comA, vec_comB):

    print(f'Zero Argument')
    print(f'Round 1: Initial message')
    com_a0, com_bm, vec_com_D = prover.zero_arg_round_1_computation(pc, params, y_had)
    
    print(f'Round 2: Verifier sends a challenge x_zero_arg:', end=" ")
    
    # Verifier sends a new challenge x
    x_zeroArg = verifier.generateChallenge(params)
    print(x_zeroArg)

    print(f'Round 3: Prover computes and sends openings')
    vec_x_m = [pow(x_zeroArg, i,q) for i in range(m+1)]
    vec_x_m_j = [pow(x_zeroArg, m-j,q) for j in range(m+1)]
    vec_x_m_k = [pow(x_zeroArg, k,q) for k in range(2*m+1)]

    vec_a, vec_b, r, s, t = prover.zero_arg_round_3_computation(params, m, x_zeroArg, vec_x_m, vec_x_m_j, vec_x_m_k)
    
    # Verification
    print(f'Round 4: Verification')
    
    verifier.zero_arg_verify(params, com_a0, com_bm, vec_com_D, vec_a, vec_b, r, s, t, vec_x_m, vec_x_m_j, vec_x_m_k, vec_comA, vec_comB, y_had )
    print(f'Zero Argument finished')   
        
'''
    comA = com(vec_a;r)
    b = prod(a_i)
    vec_a = self.vec_b_prod_arg_r1 
    comA = com_b = pc.commit(vec_b; s) 
    ComA is a commitment to rows of matrix A (d_i - z)
        Corresponds to vector: vec_b_prod_arg_r1
    b is the product of (y*i + x^i - z) for i = 1...N
    r = self.s_prod_arg_r1
'''
def singleValueProdArgument(pc, params, comA, b):
    
    print(f'Single Value Product Argument')
    print(f'Round 1: Prover`s initial message:')
    
    com_d, com_delta, com_triangle = prover.single_val_arg_round_1_computation(pc, params)

    print(f'Round 2: Verifier sends a challenge x_sing_val:', end = ' ')
    x_5_3 = verifier.generateChallenge(params)
    print(x_5_3)

    print(f'Round 3: Prover responds with openings:')
    
    vec_a_bar, vec_b_bar, r_bar, s_bar = prover.single_val_arg_round_3_computation(pc, params, x_5_3)

    # Verification
    print(f'Round 4: Verification ')
    verifier.single_value_product_arg_verify(pc, params, comA, b, com_d, com_delta, com_triangle, x_5_3, vec_a_bar, vec_b_bar, r_bar, s_bar)
    print(f'End of single value product argument')


'''4.2: Trading computation for interaction'''
def tradingCompForInter(pc, elG, m, n, params, prover, verifier, org_ciphertexts, shuffled_ciphertexts, com_A, x_init, mew):
    m_prime = m//mew 
    start = timer()
    q, _, _, p = params 
    N = m*n
    rho, vec_x = prover.multi_exponentiation_round_0_computation(params, x_init, N)

    singleCiphertext = bilinearMap(org_ciphertexts, vec_x, p)
    print(f'Product of original ciphertexts raised to vector x: {singleCiphertext}')
    
    
    print(f'Trading computation for interaction Argument:')
    print(f'Round 1: Initial Message')

    com_vec_b, vec_E = prover.multi_expo_optimization_round_1_computation(pc, elG, params, m_prime, n, shuffled_ciphertexts, rho, mew)
    
    ''''Vec_E here is a vector containing 2*mew -1 elements (where each element is a ciphertext (X,Y))'''

    print(f'Round 2: Verifier sends a challenge x_multi_optim')
    x_multi_optim = verifier.generateChallenge(params)
    print(f'Challenge sent by verifier for MEA : {x_multi_optim}')

    print(f'Round 3: Prover calculates and sends the respective openings')
    b,s, single_ciphertext_prime, new_mat_c, com_A_prime = prover.multi_expo_optimization_round_3_computation(params, m_prime, n, x_multi_optim, mew, shuffled_ciphertexts, com_A, vec_E)

    multiExpoArgumentOptimized(pc, elG, m_prime, n, params, prover, verifier, single_ciphertext_prime, new_mat_c, com_A_prime, mew )
    

    print(f'Round 4: Verification')

    verifier.multi_expo_optimization_arg_verify(pc, params, mew, singleCiphertext,  com_vec_b, vec_E, b, s, x_multi_optim)
    
    print(f'End of Multi Exponentiation Optimization Argument')
    end = timer()
    print(f'Total time taken for multi-exponentiation optimization argument: {end - start}')

'''Main argument'''
def shuffleArgument(pc, elG, params,  m, n, prover, verifier, org_ciphertexts, shuffled_ciphertexts, mew):
    start =  timer()

    print(f'Round 1: Prover sends a commitment to the permutation used')
    com_a = prover.round_1_computation(pc, params, m)

    print(f'Commitment to permutation is given by: {com_a}')
    
    print(f'Round 2: Verifier sends a challenge x:')
    x_init = verifier.generateChallenge(params)
    print(f'Challenge x_org given by: {x_init}')

    print(f'Round 3: Prover sends commitment to challenge x raised to the permutation')
    com_b = prover.round_3_computation(pc, params, m, x_init)
    print(f'Commitment to challenge raised to the permutation is given by: {com_b}')
    
    print(f'Round 4: Verifier sends two challenges')
    y_init = verifier.generateChallenge(params)
    z_init = verifier.generateChallenge(params)

    print(f'Challenge y given by: {y_init}')
    print(f'Challenge z given by: {z_init}')

    print(f'Round 5: The prover uses the product argument to show that the same permutation has been used for both commitments')
    N = m*n
    vec_z = [-z_init] * N
    vec_zero = [0] * m
    com_z = pc.generalCommit(params, vec_z, vec_zero)   # Both Prover & Verifier should be able to compute this
    print(f'Commitment to -Z given by: {com_z}')
    com_a_pow_y = []
    for i in com_a:
        com_a_pow_y.append(pow(i, y_init, p))
    com_d = [(a*b) % p for a, b in zip(com_a_pow_y, com_b)]
    print(f'Commitment to D given by: {com_d}')
    
    com_d_minus_z = prover.round_5_computation(pc, params, m, n, y_init, vec_z)

    print(f' com_ck(d - z; t ) : {com_d_minus_z}')

    # Computing the polynomial that exists in Z_q (of degree d = N)
    polyVal = 1
    for i in range(1, N+1):
        polyVal = (polyVal * ( ( ( (y_init*i) % q ) + pow(x_init, i, q) - z_init) % q)) % q
    print(f'prod(yi+x^i - z) for i = 1..16: {polyVal}')
    print(f'Product argument')
    

    
    productArgument(pc, elG, m, n, params, prover, verifier, b = polyVal, vec_comA = com_d_minus_z)

    
    # Now for Multi-Expo Argument
    multiExpoArgument(pc, elG, m, n, params, prover, verifier, org_ciphertexts, shuffled_ciphertexts, com_b, x_init )
    
    end = timer()
    totalTime = end - start 
    print(f'Shuffle Argument complete')
    print(f'Time taken by prover: {prover.prover_time()}')
    print(f'Time taken by verifier: {verifier.verifier_time()}')
    print(f'Total time taken: {totalTime} ')

    
    print(f'4.2: Multi-exponentiation with optimization')
    tradingCompForInter(pc, elG, m, n, params, prover, verifier, org_ciphertexts, shuffled_ciphertexts, com_b, x_init, mew)

pc = PedersenCommitment()
m = 4
n = 4
mew = 2       # For optimization
N = m*n  # Number of messages
prover = Prover()
verifier = Verifier()

params = pc.setup(2048, 256, N)
# The group order (q). N generators (g_i), Another generator (h), modulo (p)
q, g_i, h, p= params
g = random.choice(list(filter(lambda e: e != 1, g_i)))
x = secretsGenerator.randint(2, q-1)
print(f'g = {g}')
print(f'g**(q) mod p should be 1: {pow(g,q,p)}')
print(f'Values of generators g_i used: {g_i}')
print(f'Value of generator h used: {h}')
elG = ElGamal.construct((p, g, pow(g, x, p), x))
ElG_key_params = elG.publickey()

print(f'Prime order used: {ElG_key_params.p}')
# Generated prime satisfies p | q - 1 as well

print(f'Generator used: {ElG_key_params.g}')

print(f'Y = G^x (for private key x): {ElG_key_params.y}')


# Step 2: Convert messages to ciphertexts
'''
Each message is an element of the group G. Each ciphertext is an element of H: G x G
G itself is a subgroup of Z_p* with order q'''
messageListPowers = list_of_rand_values(N, q)

messageList = [pow(g, power, p) for power in messageListPowers]

vec_rho_org = list_of_rand_values(len(messageList), q)
org_ciphertexts = [elG._encrypt(messageList[i], vec_rho_org[i]) for i in range(len(messageList))]

print(
    f'Are all ciphertexts valid components of H (G x G)?: {groupElement(params, org_ciphertexts)}')



# Prover generates pi and rho
permuted_list, vec_rho = prover.witness_initial()

shuffled_ciphertexts = []
for k in range(N):
    shuffled_ciphertexts.append(tuple(
        (i*j) % p for i, j in zip(org_ciphertexts[permuted_list[k]-1], elG._encrypt(1, vec_rho[k]))))
print(
    f'Are shuffled ciphertexts still valid components of H: G x G?: {groupElement(params, shuffled_ciphertexts)} ')

shuffleArgument(pc, elG, params,  m, n, prover, verifier, org_ciphertexts, shuffled_ciphertexts, mew)


# To generate incorrect outputs
'''
shuffled_ciphertexts[-1] = (123, 456)       # Not a group element 
shuffleArgument(pc, elG, params,  m, n, prover, verifier, org_ciphertexts, shuffled_ciphertexts, mew)

shuffled_ciphertexts[-1] = org_ciphertexts[-1]  $ Multi-expo will fail
shuffleArgument(pc, elG, params,  m, n, prover, verifier, org_ciphertexts, shuffled_ciphertexts, mew)
'''
print('The end')
