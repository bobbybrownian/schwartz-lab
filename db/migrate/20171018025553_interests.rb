class Interests < ActiveRecord::Migration[5.1]
  def change
    create_table :interests do |t|
      t.string :description, null: false
      t.string :link, null: false
      t.integer :interest_type, null: false

      t.timestamp null: false
    end
  end
end
