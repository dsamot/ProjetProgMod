/**
 * Lire un historique depuis un fichier
 *
 * @param nbAssets nombre de sous jacents
 * @param[out] market Matrice contenant un historique de marché en sortie. Cette matrice doit
 * exister (possiblement de taille 0 x 0) avant l'appel à cette fonction
 */
void simul_market(int nbAssets, PnlMat *market)
{
    char *marketFile = "market.dat";
    PnlMat *market_from_file = pnl_mat_create_from_file(marketFile);
    pnl_mat_extract_subblock(market, market_from_file, 0, market_from_file->m, 0, nbAssets);
    pnl_mat_free(&market_from_file);
}

