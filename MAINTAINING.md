# Maintenance du dépôt source CROCO

Le moteur maintenu localement se trouve dans le dépôt séparé
`/Users/matthieucaillaud/Documents/DEV/OCEANO/croco-source`.

## Dépôts distants

- `origin` : dépôt privé GitHub `Matt29/croco` ;
- `upstream` : dépôt officiel `https://gitlab.inria.fr/croco-ocean/croco.git`.

Le dépôt part du tag officiel `v2.1.3`. Le répertoire officiel `BENCH` est
volontairement absent : les benchmarks, configurations et résultats sont
maintenus dans `AUDIT_CROCO`.

## Mettre à jour depuis CROCO officiel

Synchroniser de préférence à chaque nouvelle version stable, et avant une
évolution importante du moteur :

```bash
cd /Users/matthieucaillaud/Documents/DEV/OCEANO/croco-source
git fetch upstream --tags
git switch main
git switch -c sync/vX.Y.Z
git merge --no-ff vX.Y.Z
```

La suppression locale de `BENCH` doit être conservée lors de la fusion. Si la
version officielle modifie aussi `BENCH`, Git peut signaler des conflits
`modify/delete` ; conserver la suppression locale pour ces chemins.

Après compilation et validation des cas d'usage dans `AUDIT_CROCO`, fusionner
la branche de synchronisation dans `main`, puis publier :

```bash
git switch main
git merge --no-ff sync/vX.Y.Z
git push origin main --tags
```

## Règle de périmètre

Le dépôt source accepte les corrections et évolutions génériques du moteur :
Fortran, bibliothèques embarquées et scripts de compilation génériques. Il ne
doit pas recevoir de fichiers NetCDF, de résultats, de configurations de cas
réels, de paramètres Datarmor, ni de scripts de campagnes spécifiques.
