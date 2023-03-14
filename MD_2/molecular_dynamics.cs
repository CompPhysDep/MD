using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace molecular_dynamics
{
    class Program
    {
        public void save_system_parameters()
        {
            string st;
            FileStream fs_T = new FileStream("system_temperature.txt", FileMode.Append);
            //saving temperature
            StreamWriter sw_T = new StreamWriter(fs_T);
            st = time.ToString() + " " + temperature.ToString();
            st = st.Replace(",", ".");
            sw_T.WriteLine(st);
            sw_T.Close(); fs_T.Close();
            FileStream fs_P = new FileStream("system_pressure.txt", FileMode.Append);
            StreamWriter sw_P = new StreamWriter(fs_P);
            st = time.ToString() + " " + pressure.ToString();
            st = st.Replace(",", ".");
            sw_P.WriteLine(st);
            sw_P.Close(); fs_P.Close();
            // saving kinetic energy
            FileStream fs_Ek = new FileStream("system_kinetic_energy.txt", FileMode.Append);
            StreamWriter sw_Ek = new StreamWriter(fs_Ek);
            st = time.ToString() + " " + (kinetic_energy * N).ToString();
            st = st.Replace(",", ".");
            sw_Ek.WriteLine(st);
            sw_Ek.Close(); fs_Ek.Close();
            // potential energy
            FileStream fs_Ep = new FileStream("system_potential_energy.txt", FileMode.Append);
            StreamWriter sw_Ep = new StreamWriter(fs_Ep);
            st = time.ToString() + " " + potential_energy.ToString();
            st = st.Replace(",", ".");
            sw_Ep.WriteLine(st);
            sw_Ep.Close(); fs_Ep.Close();
            // total energy
            FileStream fs_H = new FileStream("system_total_energy.txt", FileMode.Append);
            StreamWriter sw_H = new StreamWriter(fs_H);
            st = time.ToString() + " " + H_energy.ToString();
            st = st.Replace(",", ".");
            sw_H.WriteLine(st);
            sw_H.Close(); fs_H.Close();
        }

        #region global variables
        bool move_test, update_test = true;
        float[] pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, acc_x, acc_y, acc_z;
        int N, nx, ny, nz, nb, step, ri;
        float[] e_pot;
        float H_energy;
        Random rand = new Random();
        float L_x, L_y, L_z, alat, dt, time, lambda, volume;
        int M, init_type, table_size;
        int max_list_size, list_size, J;
        string cryst_type;
        float pressure, density, virial, r_cut, rij_x, rij_y, rij_z;
        float rij2, dP, dF, ds, ds1, ds2, dr_s;//rij2 - squared radius-vector
        float[,] r_cell, table_pot;
        float[] disp_list_x, disp_list_y, disp_list_z, system_center;
        float T_ext, T0, kinetic_energy, temperature, vel_dT,
            potential_energy;/*temperature - instantaneous system temperature*/
        int[] Advance, Marker1, Marker2, List;
        #endregion
        public void load_lammps_trajectory()
        {
            string stx, sty, stz, st;
            float bx, by, bz;
            // creating vmd compatible file
            FileStream fs = new FileStream("system_config_step1.txt", FileMode.Append);
            StreamWriter wfs = new StreamWriter(fs);
            wfs.WriteLine("ITEM: TIMESTEP");
            wfs.WriteLine(1);
            wfs.WriteLine("ITEM: NUMBER OF ATOMS");
            wfs.WriteLine(N);
            wfs.WriteLine("ITEM: BOX BOUNDS pp pp pp");
            bx = L_x / 2.0f;
            by = L_y / 2.0f;
            bz = L_z / 2.0f;
            stx = (-bx).ToString() + " " + bx.ToString(); stx = stx.Replace(",", ".");
            sty = (-by).ToString() + " " + by.ToString(); sty = sty.Replace(",", ".");
            stz = (-bz).ToString() + " " + bz.ToString(); stz = stz.Replace(",", ".");
            wfs.WriteLine(stx);
            wfs.WriteLine(sty);
            wfs.WriteLine(stz);
            // record to file
            wfs.WriteLine("ITEM: ATOMS id type x y z");
            for (int i = 0; i < N; i++)
            {
                st = (i + 1).ToString() + " 0 " + pos_x[i].ToString() + " " +
                pos_y[i].ToString() + " " + pos_z[i].ToString();
                st = st.Replace(",", ".");
                wfs.WriteLine(st);
            }
            wfs.Close(); fs.Close();
        }
        /*note: if vmd throws an error\exception, increase alat parameter*/
        public void save_config()
        {
            //file for recording last step
            FileStream fss = new FileStream("system_config_last_step.txt", FileMode.Create);
            StreamWriter sww = new StreamWriter(fss);
            // record parameters of the system: Lx, Ly, Lz, N, time
            sww.WriteLine(L_x.ToString());
            sww.WriteLine(L_y.ToString());
            sww.WriteLine(L_z.ToString());
            sww.WriteLine(N.ToString());
            sww.WriteLine(time.ToString());
            for (int i = 0; i < N; i++)
            {
                // saving coordinates
                sww.WriteLine(pos_x[i].ToString() + " " + pos_y[i].ToString() + " " +
                pos_z[i].ToString());
                // saving velocities
                sww.WriteLine(vel_x[i].ToString() + " " + vel_y[i].ToString() + " " +
                vel_z[i].ToString());
                //saving accelerates
                sww.WriteLine(acc_x[i].ToString() + " " + acc_y[i].ToString() + " " + acc_z[i].ToString());
            }
            sww.Close(); fss.Close();
        }
        //boundary conditions
        public void boundary_conditions()
        {
            for (int i = 0; i < N; i++)
            {
                if (pos_x[i] >= L_x / 2.0f) pos_x[i] = pos_x[i] - L_x;
                if (pos_x[i] < -L_x / 2.0f) pos_x[i] = pos_x[i] + L_x;
                if (pos_y[i] >= L_y / 2.0f) pos_y[i] = pos_y[i] - L_y;
                if (pos_y[i] < -L_y / 2.0f) pos_y[i] = pos_y[i] + L_y;
                if (pos_z[i] >= L_z / 2.0f) pos_z[i] = pos_z[i] - L_z;
                if (pos_z[i] < -L_z / 2.0f) pos_z[i] = pos_z[i] + L_z;
            }
        }

        //searching the center of system
        public void load_system_center()
        {
            system_center = new float[3];
            system_center[0] = 0.0f;
            system_center[1] = 0.0f;
            system_center[2] = 0.0f;
            for (int i = 0; i < N; i++)
            {
                system_center[0] += pos_x[i];
                system_center[1] += pos_y[i];
                system_center[2] += pos_z[i];
            }
            system_center[0] = system_center[0] / N;
            system_center[1] = system_center[1] / N;
            system_center[2] = system_center[2] / N;
            for (int j = 0; j < N; j++)
            {
                pos_x[j] = pos_x[j] - system_center[0];
                pos_y[j] = pos_y[j] - system_center[1];
                pos_z[j] = pos_z[j] - system_center[2];
            }
        }

        public void MD()
        {
            load_initial_stat();
            load_LJ_forces();
            for (step = 1; step <= M; step++)
            {
                time = time + dt;
                system_dynamics(1);
                boundary_conditions();
                load_forces();
                system_dynamics(2);
                save_system_parameters();
                load_lammps_trajectory();

                Console.WriteLine("time = " + time + " steps/M = " +
                step.ToString() + "/" + M.ToString() + " T = " + temperature + " P = " + pressure +
                " U = " + H_energy);

            }
            save_config();
        }


        static void Main(string[] args)
        {
            Program p = new Program();
            p.MD();
            Console.ReadLine();
        }

        //creating initial conditions of the system
        public void load_initial_stat()
        {
            #region initial conditions
            nx = 7; ny = 7; nz = 7; //number of atoms along each axis
            M = 100;//steps number
            alat = 1.75f;//lattice constant
            init_type = 1;//initial type of creating system
            cryst_type = "fcc";
            dr_s = 0.3f;
            vel_dT = 0.02f;//velocity of heating/coolng, >0 is heating 
            dt = 0.005f;//timestep
            T_ext = 2.4f;//thermostat temperature
            T0 = 0.0f;//initial temperature
            #endregion
            //cell cize
            L_x = alat * nx;
            L_y = alat * ny;
            L_z = alat * nz;
            time = 0;
            if (init_type == 1) gen_order_stat(cryst_type);
            #region arrays
            vel_x = new float[N];
            vel_y = new float[N];
            vel_z = new float[N];
            acc_x = new float[N];
            acc_y = new float[N];
            acc_z = new float[N];
            e_pot = new float[N];
            Advance = new int[N];
            Marker1 = new int[N];
            Marker2 = new int[N];
            List = new int[100 * N];
            disp_list_x = new float[N];
            disp_list_y = new float[N];
            disp_list_z = new float[N];
            #endregion
        }
        //type of lattice
        public void gen_order_stat(string ctype)
        {
            if (ctype == "fcc")
            {
                // face center cubic system generation
                nb = 4;
                r_cell = new float[3, nb];
                r_cell[0, 0] = 0.0f; r_cell[1, 0] = 0.0f; r_cell[2, 0] = 0.0f;
                r_cell[0, 1] = 0.5f; r_cell[1, 1] = 0.5f; r_cell[2, 1] = 0.0f;
                r_cell[0, 2] = 0.0f; r_cell[1, 2] = 0.5f; r_cell[2, 2] = 0.5f;
                r_cell[0, 3] = 0.5f; r_cell[1, 3] = 0.0f; r_cell[2, 3] = 0.5f;
                N = 4 * nx * ny * nz;
                pos_x = new float[N];
                pos_y = new float[N];
                pos_z = new float[N];
            }
            if (ctype == "sc")
            {
                // simple cubic system generation
                nb = 1;
                r_cell = new float[3, nb];
                r_cell[0, 0] = 0.0f; r_cell[1, 0] = 0.0f; r_cell[2, 0] = 0.0f;
                N = nx * ny * nz;
                pos_x = new float[N];
                pos_y = new float[N];
                pos_z = new float[N];
            }
            int n = 0;
            for (int k = 0; k < nz; k++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int i = 0; i < nx; i++)
                    {
                        for (int L = 0; L < nb; L++)
                        {
                            n++;
                            pos_x[n - 1] = alat * (i + r_cell[0, L]);
                            pos_y[n - 1] = alat * (j + r_cell[1, L]);
                            pos_z[n - 1] = alat * (k + r_cell[2, L]);

                        }
                    }
                }
            }
        }
        //tabulated potantial of the energy and forces - Lennard-Jones case
        public void load_LJ_forces()
        {
            float sigma = 1.0f, epsilon = 1.0f, displ = 0.0001f;
            double rij = 0.5f, energy, force, energy_rcut;
            r_cut = 2.5f;
            table_size = 25000;
            table_pot = new float[table_size, 2];
            energy_rcut = 4.0 * epsilon * (Math.Pow((sigma / r_cut), 12) -
            Math.Pow((sigma / r_cut), 6));
            for (int i = 0; i < table_size; i++)
            {
                rij = rij + displ;
                energy = 4.0 * epsilon * (Math.Pow((sigma / rij), 12) -
                Math.Pow((sigma / rij), 6)) - energy_rcut;
                force = 24.0 * epsilon * (2.0 * Math.Pow((sigma / rij), 13) -
                Math.Pow((sigma / rij), 7)) / sigma;
                table_pot[i, 0] = (float)energy;
                table_pot[i, 1] = (float)force;
            }
        }


        public void system_dynamics(int stage_id)
        {
            if (stage_id == 1)
            {
                for (int i = 0; i < N; i++)
                {
                    pos_x[i] += vel_x[i] * dt + acc_x[i] * dt * dt / 2.0f;
                    pos_y[i] += vel_y[i] * dt + acc_y[i] * dt * dt / 2.0f;
                    pos_z[i] += vel_z[i] * dt + acc_z[i] * dt * dt / 2.0f;
                }
                //calculating instantaneous temperature and energy of the system
                kinetic_energy = 0;
                for (int i = 0; i < N; i++)
                {
                    kinetic_energy += vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i] + vel_z[i] * vel_z[i];
                }
                kinetic_energy = kinetic_energy / 2.0f;
                temperature = 2.0f * kinetic_energy / (3.0f * N);
                T0 = T0 + (float)(Math.Sign(T_ext - T0) * vel_dT * dt);
                //checking initial temperature is 0
                if (temperature == 0) temperature = T0;
                lambda = (float)Math.Sqrt(T0 / temperature);
                for (int i = 0; i < N; i++)
                {
                    vel_x[i] = lambda * vel_x[i] + acc_x[i] * dt / 2.0f;
                    vel_y[i] = lambda * vel_y[i] + acc_y[i] * dt / 2.0f;
                    vel_z[i] = lambda * vel_z[i] + acc_z[i] * dt / 2.0f;
                }
            }
            if (stage_id == 2)
            {
                volume = L_x * L_y * L_z;
                density = N / volume;
                pressure = density * temperature + virial / volume;
                for (int i = 0; i < N; i++)
                {
                    vel_x[i] = vel_x[i] + acc_x[i] * dt / 2.0f;
                    vel_y[i] = vel_y[i] + acc_y[i] * dt / 2.0f;
                    vel_z[i] = vel_z[i] + acc_z[i] * dt / 2.0f;
                }
                H_energy = kinetic_energy + potential_energy;
            }
        }

        public void load_forces()
        {
            for (int i = 0; i < N; i++)
            {
                acc_x[i] = 0.0f; acc_y[i] = 0.0f; acc_z[i] = 0.0f; e_pot[i] = 0.0f;
            }
            virial = 0.0f;
            if (update_test)
            {//updating Verlet
                update_list(r_cut + dr_s);
                update_test = false;
            }
            for (int i = 0; i < N - 1; i++)
            {
                for (int L = Marker1[i]; L <= Marker2[i]; L++)
                {
                    J = List[L];
                    rij_x = pos_x[i] - pos_x[J];
                    rij_y = pos_y[i] - pos_y[J];
                    rij_z = pos_z[i] - pos_z[J];
                    rij_x = rij_x - Convert.ToInt32(rij_x / L_x) * L_x;
                    rij_y = rij_y - Convert.ToInt32(rij_y / L_y) * L_y;
                    rij_z = rij_z - Convert.ToInt32(rij_z / L_z) * L_z;
                    rij2 = rij_x * rij_x + rij_y * rij_y + rij_z * rij_z;
                    if (rij2 < (r_cut * r_cut))
                    {
                        ri = (int)((Math.Pow(rij2, 0.5f) - 0.5f) / 0.0001f);
                        dP = table_pot[ri - 1, 0];
                        e_pot[i] = e_pot[i] + (dP / 2.0f);
                        e_pot[J] = e_pot[J] + (dP / 2.0f);
                        dF = table_pot[ri - 1, 1];
                        virial = virial + dF * rij2;
                        acc_x[i] = acc_x[i] + dF * rij_x;
                        acc_y[i] = acc_y[i] + dF * rij_y;
                        acc_z[i] = acc_z[i] + dF * rij_z;
                        acc_x[J] = acc_x[J] - dF * rij_x;
                        acc_y[J] = acc_y[J] - dF * rij_y;
                        acc_z[J] = acc_z[J] - dF * rij_z;
                    }
                }
            }
            virial = virial / 3.0f;
            potential_energy = 0.0f;
            for (int i = 0; i < N; i++)
            {
                potential_energy += e_pot[i];
            }
            for (int i = 0; i < N; i++)
            {
                disp_list_x[i] += vel_x[i] * dt + acc_x[i] * dt * dt / 2.0f;
                disp_list_y[i] += vel_y[i] * dt + acc_y[i] * dt * dt / 2.0f;
                disp_list_z[i] += vel_z[i] * dt + acc_z[i] * dt * dt / 2.0f;
            }
            load_move_test(dr_s);
            update_test = move_test;
        }
        public void update_list(float rang)
        {
            max_list_size = N;
            float rng;
            rng = rang * rang;
            int L = 1;
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    rij_x = pos_x[i] - pos_x[j];
                    rij_y = pos_y[i] - pos_y[j];
                    rij_z = pos_z[i] - pos_z[j];
                    rij_x = rij_x - Convert.ToInt32(rij_x / L_x) * L_x;
                    rij_y = rij_y - Convert.ToInt32(rij_y / L_y) * L_y;
                    rij_z = rij_z - Convert.ToInt32(rij_z / L_z) * L_z;
                    rij2 = rij_x * rij_x + rij_y * rij_y + rij_z * rij_z;
                    if (rij2 < rng)
                    {
                        Advance[j] = 1;
                    }
                    else
                    {
                        Advance[j] = 0;
                    }
                }
                Marker1[i] = L;
                for (int j = i + 1; j < N; j++)
                {
                    List[L] = j;
                    L = L + Advance[j];
                }
                Marker2[i] = L - 1;
            }
            list_size = L - 1;
            for (int i = 0; i < N; i++)
            {
                disp_list_x[i] = 0.0f;
                disp_list_y[i] = 0.0f;
                disp_list_z[i] = 0.0f;
            }
        }
        public void load_move_test(float skin)
        {
            ds1 = 0.0f; ds2 = 0.0f;
            for (int i = 0; i < N; i++)
            {
                ds = (float)Math.Sqrt(disp_list_x[i] * disp_list_x[i] +
                disp_list_y[i] * disp_list_y[i] + disp_list_z[i] * disp_list_z[i]);
                if (ds >= ds1)
                {
                    ds2 = ds1;
                    ds1 = ds;
                }
                else
                if (ds >= ds2)
                {
                    ds2 = ds;
                }
            }
            move_test = (ds1 + ds2 > skin);
        }
    }
}
