{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 97,
   "id": "3c3db81d",
   "metadata": {},
   "outputs": [],
   "source": [
    "import sys\n",
    "sys.path.append('../15Msun/work')\n",
    "import MESA_STELLAR as ms\n",
    "\n",
    "from subprocess import run\n",
    "from subprocess import call\n",
    "import subprocess\n",
    "import numpy as np\n",
    "import matplotlib.pyplot as plt"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 98,
   "id": "7181f188",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Updated inlist max_star_age:  150\n"
     ]
    }
   ],
   "source": [
    "# Rests inlist_star and 15Msun mod file to a 150 year old star\n",
    "ms.WriteMesaModel(ms.ReadMesaModel('../15Msun_Gen.mod'), '15ZAMtest.mod')\n",
    "ms.inlist_RESET_age('inlist_star', 55, 150 )"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 99,
   "id": "331f95e8",
   "metadata": {},
   "outputs": [],
   "source": [
    "##Global Variables\n",
    "timestep = 100\n",
    "runs = 0\n",
    "zone = 0\n",
    "flag = True\n",
    "end_star_age = 2000\n",
    "global_lnd = []\n",
    "global_lnT = []\n",
    "global_lnR = []\n",
    "global_time = []"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 100,
   "id": "4cbae730",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Starting Star Age:  1.5000000000000000D+02\n",
      "Current_Max_Age:  150.0\n",
      "Updated inlist max_star_age:  250.0\n",
      "Age from data:  250.0\n",
      "Current_Max_Age:  250.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  350.0\n",
      "MESA Starting:Run 2.000000\n",
      "Age from data:  350.0\n",
      "Current_Max_Age:  350.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  450.0\n",
      "MESA Starting:Run 3.000000\n",
      "Age from data:  450.0\n",
      "Current_Max_Age:  450.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  550.0\n",
      "MESA Starting:Run 4.000000\n",
      "Age from data:  550.0\n",
      "Current_Max_Age:  550.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  650.0\n",
      "MESA Starting:Run 5.000000\n",
      "Age from data:  650.0\n",
      "Current_Max_Age:  650.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  750.0\n",
      "MESA Starting:Run 6.000000\n",
      "Age from data:  750.0\n",
      "Current_Max_Age:  750.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  850.0\n",
      "MESA Starting:Run 7.000000\n",
      "Age from data:  850.0\n",
      "Current_Max_Age:  850.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  950.0\n",
      "MESA Starting:Run 8.000000\n",
      "Age from data:  950.0\n",
      "Current_Max_Age:  950.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  1050.0\n",
      "MESA Starting:Run 9.000000\n",
      "Age from data:  1050.0\n",
      "Current_Max_Age:  1050.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  1150.0\n",
      "MESA Starting:Run 10.000000\n",
      "Age from data:  1150.0\n",
      "Current_Max_Age:  1150.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  1250.0\n",
      "MESA Starting:Run 11.000000\n",
      "Age from data:  1250.0\n",
      "Current_Max_Age:  1250.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  1350.0\n",
      "MESA Starting:Run 12.000000\n",
      "Age from data:  1350.0\n",
      "Current_Max_Age:  1350.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  1450.0\n",
      "MESA Starting:Run 13.000000\n",
      "Age from data:  1450.0\n",
      "Current_Max_Age:  1450.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  1550.0\n",
      "MESA Starting:Run 14.000000\n",
      "Age from data:  1550.0\n",
      "Current_Max_Age:  1550.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  1650.0\n",
      "MESA Starting:Run 15.000000\n",
      "Age from data:  1650.0\n",
      "Current_Max_Age:  1650.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  1750.0\n",
      "MESA Starting:Run 16.000000\n",
      "Age from data:  1750.0\n",
      "Current_Max_Age:  1750.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  1850.0\n",
      "MESA Starting:Run 17.000000\n",
      "Age from data:  1850.0\n",
      "Current_Max_Age:  1850.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  1950.0\n",
      "MESA Starting:Run 18.000000\n",
      "Age from data:  1950.0\n",
      "Current_Max_Age:  1950.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  2050.0\n",
      "MESA Starting:Run 19.000000\n",
      "Age from data:  2050.0\n",
      "Final Star Age:  2050.0\n",
      "Current_Max_Age:  2050.0\n",
      "Star ages reaches max_age. Add timestep\n",
      "Updated inlist max_star_age:  2150.0\n",
      "MESA Starting:Run 20.000000\n",
      "Done\n"
     ]
    }
   ],
   "source": [
    "## Reading Starting Mode\n",
    "data1 = ms.ReadMesaModel(\"15ZAMtest.mod\")\n",
    "print(\"Starting Star Age: \", data1['star_age'])\n",
    "ms.inlist_curr_age('inlist_star',55)\n",
    "ms.inlist_add_age('inlist_star', 55, timestep)\n",
    "ms.RunMesa()\n",
    "runs += 1\n",
    "\n",
    "while(flag):\n",
    "    datacur = ms.ReadMesaModel(\"final_test1.mod\")\n",
    "    curstar_age = float(ms.D_toE(datacur['star_age']))\n",
    "    print(\"Age from data: \", curstar_age)\n",
    "    global_time.append(curstar_age)\n",
    "    global_lnd.append(datacur['lnd'])\n",
    "    global_lnT.append(datacur['lnT'])\n",
    "    global_lnR.append(datacur['lnR'])\n",
    "    ms.WriteMesaModel(datacur, '15ZAMtest.mod')\n",
    "\n",
    "    if(curstar_age > end_star_age):\n",
    "        print(\"Final Star Age: \", curstar_age)\n",
    "        flag = False\n",
    "    \n",
    "    if(curstar_age == ms.inlist_curr_age('inlist_star',55)):\n",
    "        print(\"Star ages reaches max_age. Add timestep\")\n",
    "        ms.inlist_add_age('inlist_star', 55, timestep)\n",
    "        print(\"MESA Starting:Run {:f}\".format(int(runs+1)))\n",
    "        ms.RunMesa()\n",
    "        runs += 1\n",
    "print(\"Done\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 115,
   "id": "ede48ec9",
   "metadata": {},
   "outputs": [],
   "source": [
    "## Data from the first Run ->> Figure out about the zones\n",
    "lnR = global_lnR[0].astype(float)\n",
    "lnT = global_lnT[0].astype(float)\n",
    "lnD = global_lnd[0].astype(float)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 116,
   "id": "648c411f",
   "metadata": {},
   "outputs": [],
   "source": [
    "## Data from Last Run\n",
    "lnRlast = global_lnR[runs-2].astype(float)\n",
    "lnTlast = global_lnT[runs-2].astype(float)\n",
    "lnDlast = global_lnd[runs-2].astype(float)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 120,
   "id": "949ab6c7",
   "metadata": {},
   "outputs": [
    {
     "ename": "FileNotFoundError",
     "evalue": "[Errno 2] No such file or directory: '../../output/1DTime_Output'",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mFileNotFoundError\u001b[0m                         Traceback (most recent call last)",
      "\u001b[0;32m/tmp/ipykernel_149377/985791986.py\u001b[0m in \u001b[0;36m<module>\u001b[0;34m\u001b[0m\n\u001b[1;32m      2\u001b[0m \u001b[0;31m# Get Data ready for output\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m      3\u001b[0m \u001b[0moutdir\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m'../../output/'\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m----> 4\u001b[0;31m \u001b[0;32mwith\u001b[0m \u001b[0mopen\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0moutdir\u001b[0m \u001b[0;34m+\u001b[0m \u001b[0;34m'1DTime_Output'\u001b[0m\u001b[0;34m)\u001b[0m \u001b[0;32mas\u001b[0m \u001b[0moutfile\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m      5\u001b[0m     \u001b[0mnp\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0msavetxt\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0moutfile\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mnp\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mcolumn_stack\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mtp\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mlnR\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mlnD\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mlnT\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m      6\u001b[0m \u001b[0mprint\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m\"Done\"\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;31mFileNotFoundError\u001b[0m: [Errno 2] No such file or directory: '../../output/1DTime_Output'"
     ]
    }
   ],
   "source": [
    "########################################################################################################3\n",
    "# Get Data ready for output\n",
    "outdir = '../../output/'\n",
    "with open(outdir + '1DTime_Output') as outfile:\n",
    "    np.savetxt(outfile, np.column_stack([tp, lnR, lnD, lnT]))\n",
    "print(\"Done\")"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.8"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
