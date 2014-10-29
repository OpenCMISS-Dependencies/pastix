# Nom arbitraire du travail LoadLeveler
# @ job_name = _NAME_
# Fichier de sortie standard du travail
# @ output = _FILE_OUT_
# Fichier de sortie d'erreur du travail
# @ error =  _FILE_ERR_
# Type de travail
# @ job_type = parallel
# pour recevoir un mail en cas de time limit exceeded (ou autre pb.)
# @ notification = complete
# Nombre de processus MPI demandes
# @ total_tasks = _NBPROC_
# Nombre de taches OpenMP/pthreads par processus MPI
# @ parallel_threads = _NBTHREAD_
# Memoire max. utilisee par processus
# @ data_limit = 102.4Gb
# Temps du job hh:mm:ss
# @ wall_clock_limit = _TIME_
# User to notify
# @ notify_user = _USER_
# @ queue

# Pour avoir l'echo des commandes
#set -x

# Repertoire temporaire de travail
cd _PATH_

# La memoire STACK max. (defaut 4Mo) utilisee (ici 64 Mo) par
# les variables privees de chaque thread.
export XLSMPOPTS=stack=65536000

# Execution du programme parallele mixte
_EXECCMD_  _EXECUTABLE_  _DRIVER_
