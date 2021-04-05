mi <- readr::read_csv(file = paste0(out_dir, "mi.csv"), col_names = F) %>% as.matrix()
mt <- readr::read_csv(file = paste0(out_dir, "mt.csv"), col_names = F) %>% as.matrix()
mseg <- readr::read_csv(file = paste0(out_dir, "mseg.csv"), col_names = F) %>% as.matrix()
mh <- readr::read_csv(file = paste0(out_dir, "mh.csv"), col_names = F) %>% as.matrix()

sam = NULL
vv = 1
while (vv <= num_chain) {
sam <- rbind(sam, mylist[[vv]])
vv = vv + 1
}

## distance calculation can be fully vectorized (and storing z), but I don't have time for this.
z <- sam[, stringr::str_which(cname, paste0("^z\\.[0-1]\\."))] %>% as.matrix()
w <- sam[, stringr::str_which(cname, paste0("^w\\.[0-1]\\."))] %>% as.matrix()
## beta <- sam[, stringr::str_which(cname, paste0("^beta\\.", item, "\\."))]

## correction for single z or w
if (ncol(z) < 2 * 2 * N) z <- cbind(z, z)
if (ncol(w) < 2 * 2 * I) w <- cbind(w, w)

theta <- sam[, stringr::str_which(cname, paste0("^theta\\."))] %>% as.matrix() ## (theta.k.0 theta.k.1)
gamma <- sam[, stringr::str_which(cname, paste0("gamma\\."))] %>% as.matrix()

lambda <- sam[, stringr::str_which(cname, paste0("^lambda\\.[0-1]\\."))] %>% as.matrix()

param <- list(I = as.integer(I), N = as.integer(N), G = as.integer(G), seg = as.matrix(mseg), H = as.matrix(mh), len = as.integer(diff(sj)), sj = as.numeric(sj), Y = mi)
