## parametrically generating counterbalanced presentation lists
n_item <- NULL

ivs <- list(A = c("A1", "A2"), # bsbi
            B = c("B1", "B2"),
            C = c("C1", "C2"), # bswi
            D = c("D1", "D2"),
            E = c("E1", "E2"), # wsbi
            F = c("F1", "F2"),
            G = c("G1", "G2"), #wswi
            H = c("H1", "H2"))

ivs_3 <- list(I = paste0("A", 1:3),
              J = paste0("B", 1:3),
              K = paste0("C", 1:3))

subj_between <- c("A", "B", "C", "D")
item_between <- c("A", "B", "E", "F")

ff <- pres_lists(ivs,
                 subj_between = c("A", "B", "C", "D"),
                 item_between = c("A", "B", "E", "F"))

ff <- pres_lists(ivs[1:4],
                 subj_between = "A",
                 item_between = c("A", "B"))

ff <- pres_lists(ivs[1:4],
                 item_between = c("A", "B"))

ff <- pres_lists(ivs[1:4],
                 item_between = c("A"))

ff <- pres_lists(ivs[1:4])

ff <- pres_lists(ivs[1:4],
                 subj_between = "A",
                 n_item = 24)

ff <- pres_lists(ivs[1:3],
                 subj_between = names(ivs[1:3]),
                 item_between = names(ivs[1:3]),
                 n_item = 40)

pres_lists(ivs[1], n_item = 8)

pres_lists(ivs[1],
           subj_between = names(ivs)[1],
           item_between = names(ivs)[1],
           n_item = 8)

pres_lists(ivs_3,
           subj_between = names(ivs_3),
           n_item = 10)
