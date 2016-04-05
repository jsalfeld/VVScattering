from MitAna.TreeMod.bambu import mithep

leptonExampleMod = mithep.LeptonExampleMod(
    VertexName = 'GoodVertexes',
    MuonName = 'FiducialMuons',
    ElectronName = 'FiducialElectrons',
    MuonIdName = 'BaselineMuonId',
    ElectronIdName = 'FiducialElectronId'
)
